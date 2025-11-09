#!/bin/bash
# filepath: d:\MD Sims\SME_SE\scripts\structure_gen.sh
set -euo pipefail

# ============================================
# CONFIGURATION - Edit these variables
# ============================================
PROJECT_ROOT="/home/rimuru/workspace"  # Windows path (Git Bash compatible)
SEED=42

SIZE=20           # Diameter in nm
SHAPE="blob"      # sphere, blob, rough, faceted, or ellipsoid
NI_PERCENT=50.0   # Ni percentage (45-55 typical for NiTi)

# CPU CORES (NEW: -1 = use all available cores)
CORES=-1          # Number of CPU cores (-1 = auto-detect all cores)

# GRAIN STRUCTURE OPTIONS (choose ONE):
# Option 1: Specify target grain size (RECOMMENDED for physical realism)
GRAIN_SIZE=5.0    # Target grain diameter in nm (0 = single crystal)

# Option 2: Specify exact number of grains (advanced control)
# GRAINS=64       # Uncomment to use instead of GRAIN_SIZE

# DEFECTS (optional, realistic for EDM synthesis)
VACANCIES=0.01           # Bulk vacancy concentration (1%)
SURFACE_VACANCIES=0.03   # Surface vacancy concentration (3%)
ANTISITES=0.0            # Antisite defect concentration (0%)

# ATOMSK INTEGRATION (optional, requires Atomsk installed)
USE_ATOMSK=false  # true = use Atomsk for grain boundaries, false = use Voronoi

# ============================================
# AUTO-GENERATED PARAMETERS (do not edit)
# ============================================

# Determine structure type
if (( $(python3 -c "print(int(${GRAIN_SIZE} > 0))") )); then
    TYPE="nanocrystalline_${GRAIN_SIZE}nm_grains"
elif [[ -n "${GRAINS:-}" ]] && [[ ${GRAINS} -gt 1 ]]; then
    TYPE="polycrystalline_${GRAINS}grains"
else
    TYPE="single_crystal"
fi

# Output filename
NI_INT=$(python3 -c "print(int(${NI_PERCENT}))")
OUTPUT="${PROJECT_ROOT}/data/structures/niti_${SIZE}nm_Ni${NI_INT}_${TYPE}.data"
OUTPUT_XYZ="${OUTPUT%.data}.xyz"

# Create output directory
mkdir -p "${PROJECT_ROOT}/data/structures"

echo "=========================================="
echo "GENERATING NiTi NANOPARTICLE"
echo "=========================================="
echo "Size:          ${SIZE} nm"
echo "Shape:         ${SHAPE}"
echo "Composition:   Ni${NI_INT}Ti$(python3 -c "print(100-${NI_INT})")"
echo "Structure:     ${TYPE}"
echo "Seed:          ${SEED}"
echo "Seed:          ${SEED}"
if [[ ${CORES} -eq -1 ]]; then
    DETECTED_CORES=$(python3 -c "import multiprocessing; print(multiprocessing.cpu_count())")
    echo "CPU Cores:     ${DETECTED_CORES} (all available)"
else
    echo "CPU Cores:     ${CORES}"
fi
echo "=========================================="

# Build command
CMD="python3 \"${PROJECT_ROOT}/src/ase_np_gen.py\" \
    --diameter ${SIZE} \
    --shape ${SHAPE} \
    --ni-percent ${NI_PERCENT} \
    --output \"${OUTPUT}\" \
    --seed ${SEED} \
    --cores ${CORES} \
    --xyz"

# Add grain structure (use grain-size if specified, otherwise use grains)
if (( $(python3 -c "print(int(${GRAIN_SIZE:-0} > 0))") )); then
    CMD="${CMD} --grain-size ${GRAIN_SIZE}"
    echo "Grain size:    ${GRAIN_SIZE} nm (auto-calculated grain count)"
elif [[ -n "${GRAINS:-}" ]] && [[ ${GRAINS} -gt 1 ]]; then
    CMD="${CMD} --polycrystalline ${GRAINS}"
    echo "Grains:        ${GRAINS} (exact count)"
else
    echo "Grains:        Single crystal"
fi

# Add defects
if (( $(python3 -c "print(int(${VACANCIES} > 0))") )); then
    CMD="${CMD} --vacancies ${VACANCIES}"
    echo "Vacancies:     ${VACANCIES} (bulk)"
fi

if (( $(python3 -c "print(int(${SURFACE_VACANCIES} > 0))") )); then
    CMD="${CMD} --surface-vacancies ${SURFACE_VACANCIES}"
    echo "Surf. vac.:    ${SURFACE_VACANCIES}"
fi

if (( $(python3 -c "print(int(${ANTISITES} > 0))") )); then
    CMD="${CMD} --antisites ${ANTISITES}"
    echo "Antisites:     ${ANTISITES}"
fi

# Add Atomsk flag if enabled
if [[ "${USE_ATOMSK}" == "true" ]]; then
    CMD="${CMD} --use-atomsk"
    echo "Method:        Atomsk (peer-reviewed GB generation)"
else
    echo "Method:        Voronoi tessellation"
fi

echo "=========================================="
echo "Executing: ${CMD}"
echo "=========================================="

# Execute command
eval ${CMD}

# Check if files were created
if [[ -f "${OUTPUT}" ]]; then
    echo ""
    echo "=========================================="
    echo "✓ SUCCESS!"
    echo "=========================================="
    echo "LAMMPS data:   ${OUTPUT}"
    echo "XYZ file:      ${OUTPUT_XYZ}"
    
    # Print file sizes
    LAMMPS_SIZE=$(stat -c%s "${OUTPUT}" 2>/dev/null || stat -f%z "${OUTPUT}" 2>/dev/null || echo "unknown")
    XYZ_SIZE=$(stat -c%s "${OUTPUT_XYZ}" 2>/dev/null || stat -f%z "${OUTPUT_XYZ}" 2>/dev/null || echo "unknown")
    
    echo "File sizes:"
    echo "  LAMMPS: $(python3 -c "print(f'{${LAMMPS_SIZE}/1024:.1f} KB')" 2>/dev/null || echo "${LAMMPS_SIZE} bytes")"
    echo "  XYZ:    $(python3 -c "print(f'{${XYZ_SIZE}/1024:.1f} KB')" 2>/dev/null || echo "${XYZ_SIZE} bytes")"
    echo "=========================================="
    
    # Print next steps
    echo ""
    echo "NEXT STEPS:"
    echo "1. Visualize: ovito ${OUTPUT_XYZ}"
    echo "2. Minimize:  bash scripts/minimize_structure.sh"
    echo "3. Analyze:   Check grain boundaries and defects in OVITO"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "❌ ERROR: Output file not created!"
    echo "=========================================="
    exit 1
fi