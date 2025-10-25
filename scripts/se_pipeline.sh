#!/bin/bash
# filepath: /home/rimuru/workspace/scripts/se_pipeline.sh
# MODIFIED: Runs jobs sequentially for optimal performance on 8-core VM.
set -euo pipefail

# ============================================
# CONFIGURATION
# ============================================
DATA_FILES=(
    "/home/rimuru/workspace/data/structures/niti_2nm_Ni50_amorphous.data"
    # "data/structures/niti_10nm_Ni50_amorphous.data"
    # "data/structures/niti_20nm_Ni50_crystalline_1grain.data"
    # "data/structures/niti_20nm_Ni50_amorphous.data"
    # "data/structures/niti_30nm_Ni50_crystalline_1grain.data"
    # "data/structures/niti_30nm_Ni50_amorphous.data"
)

PROJECT_ROOT="/home/rimuru/workspace"
LAMMPS_BIN="/usr/local/bin/lmp"
SEED=42

# FIXED: Define potential paths as variables
POT_LIB="/home/rimuru/workspace/potentials/library.meam"
POT_FILE="/home/rimuru/workspace/potentials/NiTi.meam"

TOTAL_CORES=$(nproc)
TOTAL_MEM_GB=$(awk '/MemTotal/ {printf "%.0f\n", $2/1024/1024}' /proc/meminfo)

# ============================================
# Pre-flight Checks
# ============================================
if [[ ! -x "${LAMMPS_BIN}" ]]; then
    echo "ERROR: LAMMPS binary not found at ${LAMMPS_BIN}" >&2
    echo "Please install LAMMPS or update LAMMPS_BIN path" >&2
    exit 1
fi

# Check for potential files
if [[ ! -f "${POT_LIB}" || ! -f "${POT_FILE}" ]]; then
    echo "ERROR: Potential file not found." >&2
    echo "Checked paths:" >&2
    echo "  ${POT_LIB}" >&2
    echo "  ${POT_FILE}" >&2
    exit 1
fi

export OMP_NUM_THREADS=1

# ============================================
# Setup
# ============================================
cd "${PROJECT_ROOT}"
mkdir -p logs

echo "=========================================="
echo "SE Pipeline: Minimize → Mechanical Loading"
echo "(Running SERIALLY)"
echo "=========================================="
echo "Hardware detected:"
echo "   • CPU cores: ${TOTAL_CORES}"
echo "   • Total RAM: ${TOTAL_MEM_GB} GB"
echo "   • Files to process: ${#DATA_FILES[@]}"
echo "=========================================="

# ============================================
# Helper Functions
# ============================================
get_base_name() {
    basename "$1" .data
}

get_atom_count() {
    local datafile=$1
    [[ ! -f "${datafile}" ]] && echo "0" && return
    local count
    count=$(grep -m1 "atoms" "${datafile}" 2>/dev/null | awk '{print $1}')
    echo "${count:-0}"
}

estimate_memory() {
    local atoms=$1
    # MEAM: ~18 KB/atom → convert to GB: ÷1024 (→MB) ÷1024 (→GB) + 2GB overhead
    # Lowered minimum from 3 to 1 for testing on small VMs
    python3 -c "print(max(1, int(${atoms} * 18 / 1024 / 1024) + 1))"
}

assign_ranks() {
    local atoms=$1
    # Do not assign more ranks than we have cores
    local max_ranks=${TOTAL_CORES}
    
    if [[ ${atoms} -lt 10000 ]]; then
        echo "2"
    elif [[ ${atoms} -lt 100000 ]]; then
        echo "4"
    elif [[ ${atoms} -lt 500000 ]]; then
        echo "8"
    else
        echo "${max_ranks}"
    fi
}

# ============================================
# STEP 1: Analyze structures
# ============================================
echo ""
echo "[1/3] Analyzing structures..."

declare -A ATOM_COUNTS
declare -A MEM_ESTIMATES
declare -A MPI_RANKS

for datafile in "${DATA_FILES[@]}"; do
    if [[ ! -f "${datafile}" ]]; then
        echo "WARNING: File not found: ${datafile}" >&2
        continue
    fi
    
    structure_name=$(get_base_name "${datafile}")
    atoms=$(get_atom_count "${datafile}")
    
    if [[ ${atoms} -eq 0 ]]; then
        echo "WARNING: Could not read atom count from ${datafile}" >&2
        continue
    fi
    
    ranks=$(assign_ranks "${atoms}")
    mem=$(estimate_memory "${atoms}")
    
    ATOM_COUNTS["${structure_name}"]=${atoms}
    MEM_ESTIMATES["${structure_name}"]=${mem}
    MPI_RANKS["${structure_name}"]=${ranks}
    
    echo "   • ${structure_name}"
    echo "       Atoms: ${atoms} | RAM: ~${mem} GB | MPI ranks: ${ranks}"
done

# ============================================
# STEP 2: Minimize structures
# ============================================
echo ""
echo "[2/3] Minimizing structures (sequentially)..."

for datafile in "${DATA_FILES[@]}"; do
    [[ ! -f "${datafile}" ]] && continue
    
    structure_name=$(get_base_name "${datafile}")
    outdir="data/raw_output/${structure_name}/minimize"
    marker="${outdir}/.minimize_complete"
    output_check="${outdir}/minimized.data"
    
    if [[ -f "${marker}" ]]; then
        echo "   ✓ ${structure_name} (already done)"
        continue
    fi
    
    mkdir -p "${outdir}"
    
    mem="${MEM_ESTIMATES[$structure_name]}"
    ranks="${MPI_RANKS[$structure_name]}"
    
    echo "   → Running ${structure_name} (${ranks} ranks, ~${mem} GB)..."
    
    logfile="logs/min_${structure_name}.log"
    
    # Run the command in the foreground (no '&')
    # FIXED: Added variables for potential files
    mpirun --allow-run-as-root -n "${ranks}" "${LAMMPS_BIN}" \
        -var structure_file "${datafile}" \
        -var output_dir "${outdir}" \
        -var pot_lib "${POT_LIB}" \
        -var pot_file "${POT_FILE}" \
        -in src/lmps/minimize.lmp \
        2>&1 | tee "${logfile}"
    
    # Check if the job was successful
    if [[ -f "${output_check}" ]]; then
        touch "${marker}"
        echo "     ✓ ${structure_name} complete."
    else
        echo "     ✗ ${structure_name} FAILED (check ${logfile})" >&2
        exit 1
    fi
done

echo "   ✓ All minimizations complete"

# ============================================
# STEP 3: SE mechanical loading
# ============================================
echo ""
echo "[3/3] Running SE mechanical loading (sequentially)..."

for datafile in "${DATA_FILES[@]}"; do
    [[ ! -f "${datafile}" ]] && continue
    
    structure_name=$(get_base_name "${datafile}")
    indata="data/raw_output/${structure_name}/minimize/minimized.data"
    outdir="data/raw_output/${structure_name}/se"
    marker="${outdir}/.se_complete"
    # Use your robust check: final.data or stress_strain.txt
    output_check_1="${outdir}/final.data"
    output_check_2="${outdir}/stress_strain.txt"
    
    if [[ -f "${marker}" ]]; then
        echo "   ✓ ${structure_name} (already done)"
        continue
    fi
    
    if [[ ! -f "${indata}" ]]; then
        echo "   ✗ ${structure_name} - minimized data missing, skipping" >&2
        continue
    fi
    
    mkdir -p "${outdir}"
    
    mem="${MEM_ESTIMATES[$structure_name]}"
    ranks="${MPI_RANKS[$structure_name]}"

    echo "   → Running ${structure_name} (${ranks} ranks, ~${mem} GB)..."
    
    logfile="logs/se_${structure_name}.log"
    
    # Run the command in the foreground (no '&')
    # FIXED: Added variables for potential files
    mpirun --allow-run-as-root -n "${ranks}" "${LAMMPS_BIN}" \
        -var input_data "${indata}" \
        -var output_dir "${outdir}" \
        -var seed ${SEED} \
        -var pot_lib "${POT_LIB}" \
        -var pot_file "${POT_FILE}" \
        -in src/lmps/se_mechanical_load.lmp \
        2>&1 | tee "${logfile}"
        
    # Check if the job was successful
    if [[ -f "${output_check_1}" || -f "${output_check_2}" ]]; then
        touch "${marker}"
        echo "     ✓ ${structure_name} complete."
    else
        echo "     ✗ ${structure_name} FAILED (check ${logfile})" >&2
        exit 1
    fi
done

echo "   ✓ SE mechanical loading complete"

# ============================================
# Summary
# ============================================
echo ""
echo "=========================================="
echo "✓ SE PIPELINE COMPLETE"
echo "=========================================="
echo "Results:"
for datafile in "${DATA_FILES[@]}"; do
    [[ ! -f "${datafile}" ]] && continue
    structure_name=$(get_base_name "${datafile}")
    if [[ -f "data/raw_output/${structure_name}/se/.se_complete" ]]; then
        echo "   ✓ ${structure_name}"
    else
        echo "   ✗ ${structure_name} (INCOMPLETE)"
    fi
done
echo ""
echo "Output: data/raw_output/${structure_name}/{minimize,se}/"
echo "Logs: logs/"
echo "=========================================="