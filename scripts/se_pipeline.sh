#!/bin/bash
set -euo pipefail

PROJECT_ROOT="/workspace"
LAMMPS_BIN="/opt/deepmd-kit/install/bin/lmp"
CPU_CORES=8
SEED=42

cd "${PROJECT_ROOT}"

mkdir -p logs

echo "=========================================="
echo "SE Pipeline: Minimize → Mechanical Loading"
echo "8 cores | 6 structures"
echo "=========================================="

# ============================================================
# STEP 1: Generate structures (if needed)
# ============================================================
echo ""
echo "[1/3] Checking structures..."

if [[ ! -f "data/structures/niti_30nm_amor.data" ]]; then
    echo "  → Structures missing, running generation..."
    bash scripts/01_generate_structures.sh
else
    echo "  ✓ All structures exist"
fi

# ============================================================
# STEP 2: Minimize (energy relaxation)
# ============================================================
echo ""
echo "[2/3] Minimizing structures..."

declare -a PIDS

run_min() {
    SIZE=$1
    TYPE=$2
    MPI_RANKS=$3
    
    STRUCT="data/structures/niti_${SIZE}nm_${TYPE}.data"
    OUTDIR="data/raw_output/${SIZE}nm_${TYPE}/minimize"
    MARKER="${OUTDIR}/.minimize_complete"
    
    # Skip if already minimized
    if [[ -f "${MARKER}" ]]; then
        echo "  ✓ ${SIZE}nm ${TYPE} (already minimized, skipping)"
        return 0
    fi
    
    mkdir -p "${OUTDIR}"
    echo "  → ${SIZE}nm ${TYPE} (${MPI_RANKS} MPI ranks)"
    
    if [[ ! -f "${STRUCT}" ]]; then
        echo "ERROR: Structure file not found: ${STRUCT}" >&2
        return 1
    fi
    
    export OMP_NUM_THREADS=1
    
    mpirun --allow-run-as-root -n "${MPI_RANKS}" "${LAMMPS_BIN}" \
        -var structure_file "${STRUCT}" \
        -var output_dir "${OUTDIR}" \
        -in src/lammps/templates/minimize.lmp \
        > "logs/min_${SIZE}nm_${TYPE}.log" 2>&1 &
    
    local pid=$!
    PIDS+=($pid)
    
    # Create marker file on successful completion
    (wait ${pid} && touch "${MARKER}") &
}

wait_batch() {
    local failed=0
    for pid in "${PIDS[@]}"; do
        if ! wait "${pid}"; then
            echo "ERROR: Job with PID ${pid} failed. Check logs/" >&2
            failed=1
        fi
    done
    PIDS=()
    
    if [[ ${failed} -eq 1 ]]; then
        exit 1
    fi
}

# Batch 1: 10nm (4 cores total)
run_min 10 sc 2
run_min 10 amor 2
wait_batch

# Batch 2: 20nm (8 cores total)
run_min 20 sc 4
run_min 20 amor 4
wait_batch

# Batch 3: 30nm sc (8 cores)
run_min 30 sc 8
wait_batch

# Batch 4: 30nm amor (8 cores)
run_min 30 amor 8
wait_batch

echo "  ✓ All minimizations complete"

# ============================================================
# STEP 3: SE mechanical loading
# ============================================================
echo ""
echo "[3/3] Running SE mechanical loading..."

run_se() {
    SIZE=$1
    TYPE=$2
    MPI_RANKS=$3
    
    INDATA="data/raw_output/${SIZE}nm_${TYPE}/minimize/minimized.data"
    OUTDIR="data/raw_output/${SIZE}nm_${TYPE}/se"
    MARKER="${OUTDIR}/.se_complete"
    
    # Skip if already done
    if [[ -f "${MARKER}" ]]; then
        echo "  ✓ ${SIZE}nm ${TYPE} (already complete, skipping)"
        return 0
    fi
    
    mkdir -p "${OUTDIR}"
    echo "  → ${SIZE}nm ${TYPE} (${MPI_RANKS} MPI ranks)"
    
    if [[ ! -f "${INDATA}" ]]; then
        echo "ERROR: Minimized data not found: ${INDATA}" >&2
        return 1
    fi
    
    export OMP_NUM_THREADS=1
    
    mpirun --allow-run-as-root -n "${MPI_RANKS}" "${LAMMPS_BIN}" \
        -var input_data "${INDATA}" \
        -var output_dir "${OUTDIR}" \
        -var seed ${SEED} \
        -in src/lammps/templates/se_mechanical_load.lmp \
        > "logs/se_${SIZE}nm_${TYPE}.log" 2>&1 &
    
    local pid=$!
    PIDS+=($pid)
    
    # Create marker file on successful completion
    (wait ${pid} && touch "${MARKER}") &
}

# Batch 1: 10nm (4 cores total)
run_se 10 sc 2
run_se 10 amor 2
wait_batch

# Batch 2: 20nm (8 cores total)
run_se 20 sc 4
run_se 20 amor 4
wait_batch

# Batch 3: 30nm sc (8 cores)
run_se 30 sc 8
wait_batch

# Batch 4: 30nm amor (8 cores)
run_se 30 amor 8
wait_batch

echo "  ✓ SE mechanical loading complete"

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "✓ SE PIPELINE COMPLETE"
echo "=========================================="
echo "Outputs:"
echo "  • Minimized:  data/raw_output/*/minimize/"
echo "  • SE results: data/raw_output/*/se/"
echo "  • Logs:       logs/"
echo ""
echo "Next: Analyze SE stress-strain curves and superelastic recovery"