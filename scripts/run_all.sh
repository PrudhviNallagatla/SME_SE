#!/bin/bash
set -euo pipefail

PROJECT_ROOT="/workspace"
LAMMPS_BIN="/opt/deepmd-kit/install/bin/lmp"
CPU_CORES=8
SEED=42

cd "${PROJECT_ROOT}"

mkdir -p data/structures logs

echo "=========================================="
echo "NiTi Nanoparticle Pipeline - Full Run"
echo "8 cores | 6 structures (10/20/30 nm × sc/amor)"
echo "=========================================="

# ============================================================
# STEP 1: Generate structures
# ============================================================
echo ""
echo "[1/3] Generating NiTi structures..."

for SIZE in 10 20 30; do
    echo "  → ${SIZE}nm (single-crystal + amorphous)"
    
    python3 src/niti_np_gen.py \
        --diameter ${SIZE} --shape blob --ni-percent 50.0 \
        --output data/structures/niti_${SIZE}nm_sc.data \
        --lammps --seed ${SEED}

    python3 src/niti_np_gen.py \
        --diameter ${SIZE} --shape blob --amorphous 0.3 --ni-percent 50.0 \
        --output data/structures/niti_${SIZE}nm_amor.data \
        --lammps --seed ${SEED}
done

echo "✓ Generated 6 structures"

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
    
    PIDS+=($!)
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

echo "✓ Minimized all structures"

# ============================================================
# STEP 3: SME thermal cycling
# ============================================================
echo ""
echo "[3/3] Running SME thermal cycling..."

run_sme() {
    SIZE=$1
    TYPE=$2
    MPI_RANKS=$3
    
    INDATA="data/raw_output/${SIZE}nm_${TYPE}/minimize/minimized.data"
    OUTDIR="data/raw_output/${SIZE}nm_${TYPE}/sme"
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
        -in src/lammps/templates/sme_thermal_cycle.lmp \
        > "logs/sme_${SIZE}nm_${TYPE}.log" 2>&1 &
    
    PIDS+=($!)
}

# Batch 1: 10nm (4 cores total)
run_sme 10 sc 2
run_sme 10 amor 2
wait_batch

# Batch 2: 20nm (8 cores total)
run_sme 20 sc 4
run_sme 20 amor 4
wait_batch

# Batch 3: 30nm sc (8 cores)
run_sme 30 sc 8
wait_batch

# Batch 4: 30nm amor (8 cores)
run_sme 30 amor 8
wait_batch

echo "✓ SME thermal cycling complete"

# ============================================================
# Summary
# ============================================================
echo ""
echo "=========================================="
echo "✓ PIPELINE COMPLETE"
echo "=========================================="
echo "Outputs:"
echo "  • Structures:  data/structures/"
echo "  • Minimized:   data/raw_output/*/minimize/"
echo "  • SME results: data/raw_output/*/sme/"
echo "  • Logs:        logs/"
echo ""
echo "Next: Analyze SME results (shape recovery, hysteresis, etc.)"