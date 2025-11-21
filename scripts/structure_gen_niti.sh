#!/bin/bash
set -euo pipefail

# ============================================
# CONFIGURATION - SCIENTIFICALLY VALIDATED
# ============================================
PROJECT_ROOT="/home/rimuru/workspace"
SEED=42

SIZE=20
SHAPE="blob"
NI_PERCENT=50.0

GRAINS=0
AMORPHOUS=0.0

VACANCIES=0.01
SURFACE_VACANCIES=0.03
ANTISITES=0.005

# ============================================
# PYTHON VALIDATION (AUTO-CHECK)
# ============================================
echo "=========================================="
echo "VALIDATING SCIENTIFIC TOOLS"
echo "=========================================="

# Validate Python dependencies
if ! python3 -c "import numpy" 2>/dev/null; then
    echo "❌ ERROR: Missing Python dependencies"
    echo "   Run: pip install numpy"
    exit 1
fi
echo "✅ Python dependencies validated"

echo "=========================================="

# ============================================
# STRUCTURE GENERATION - IMPROVED VALIDATION
# ============================================

# Validate AMORPHOUS range
if (( $(python3 -c "print(int(${AMORPHOUS} < 0 or ${AMORPHOUS} > 1))") )); then
    echo "❌ ERROR: AMORPHOUS must be between 0.0 and 1.0"
    echo "   Current value: ${AMORPHOUS}"
    exit 1
fi

# Validate structure parameters (IMPROVED CHECK)
if [ "${GRAINS}" -gt 0 ] && (( $(python3 -c "print(int(${AMORPHOUS} > 0))") )); then
    echo "❌ ERROR: Cannot use both GRAINS and AMORPHOUS!"
    echo "   GRAINS=${GRAINS} (polycrystalline)"
    echo "   AMORPHOUS=${AMORPHOUS} (disordered)"
    echo ""
    echo "These are MUTUALLY EXCLUSIVE structure types."
    echo "Choose ONLY ONE:"
    echo "  - Set GRAINS=0 and AMORPHOUS=0.0 for single crystal"
    echo "  - Set GRAINS=N (N>0) and AMORPHOUS=0.0 for polycrystalline"
    echo "  - Set GRAINS=0 and AMORPHOUS=0.3 for amorphous (EDM)"
    exit 1
fi

# Determine structure type
if (( $(python3 -c "print(int(${AMORPHOUS} > 0))") )); then
    TYPE="amorphous_disorder${AMORPHOUS}"
    STRUCTURE_TYPE="amorphous"
elif (( $(python3 -c "print(int(${GRAINS} > 0))") )); then
    TYPE="polycrystalline_${GRAINS}_grains"
    STRUCTURE_TYPE="polycrystalline"
else
    TYPE="single_crystal"
    STRUCTURE_TYPE="single_crystal"
fi

# Output filename
NI_INT=$(python3 -c "print(int(${NI_PERCENT}))")
OUTPUT="${PROJECT_ROOT}/data/structures/niti_${SIZE}nm_Ni${NI_INT}_${TYPE}.data"
OUTPUT_XYZ="${OUTPUT%.data}.xyz"

mkdir -p "${PROJECT_ROOT}/data/structures"

echo "=========================================="
echo "GENERATING NiTi NANOPARTICLE"
echo "=========================================="
echo "Size:          ${SIZE} nm"
echo "Shape:         ${SHAPE}"
echo "Composition:   Ni${NI_INT}Ti$(python3 -c "print(100-${NI_INT})")"
echo "Structure:     ${STRUCTURE_TYPE}"

if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
    echo "Grains:        ${GRAINS}"
elif [[ "${STRUCTURE_TYPE}" == "amorphous" ]]; then
    echo "Disorder:      ${AMORPHOUS} (rapid quench simulation)"
fi

echo "Seed:          ${SEED}"
echo "=========================================="

# Build command
CMD="python3 \"${PROJECT_ROOT}/src/niti_np_gen.py\" \
    --diameter ${SIZE} \
    --shape ${SHAPE} \
    --ni-percent ${NI_PERCENT} \
    --output \"${OUTPUT_XYZ}\" \
    --seed ${SEED} \
    --lammps"

# ✅ Add structure type (ONLY ONE will be active)
if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
    CMD="${CMD} --grains ${GRAINS}"
elif [[ "${STRUCTURE_TYPE}" == "amorphous" ]]; then
    CMD="${CMD} --amorphous ${AMORPHOUS}"
fi

# Add defects
CMD="${CMD} --vacancies ${VACANCIES}"
CMD="${CMD} --surface-vac ${SURFACE_VACANCIES}"
CMD="${CMD} --antisites ${ANTISITES}"

echo "Executing: ${CMD}"
echo "=========================================="

# Execute
eval ${CMD}

# Final validation
if [[ -f "${OUTPUT}" ]]; then
    echo ""
    echo "=========================================="
    echo "✅ GENERATION SUCCESSFUL"
    echo "=========================================="
    echo "LAMMPS data:   ${OUTPUT}"
    echo "XYZ file:      ${OUTPUT_XYZ}"
    
    LAMMPS_SIZE=$(stat -c%s "${OUTPUT}" 2>/dev/null || stat -f%z "${OUTPUT}" 2>/dev/null || echo "0")
    XYZ_SIZE=$(stat -c%s "${OUTPUT_XYZ}" 2>/dev/null || stat -f%z "${OUTPUT_XYZ}" 2>/dev/null || echo "0")
    
    echo ""
    echo "File sizes:"
    echo "  LAMMPS: $(python3 -c "print(f'{${LAMMPS_SIZE}/1024:.1f} KB')")"
    echo "  XYZ:    $(python3 -c "print(f'{${XYZ_SIZE}/1024:.1f} KB')")"
    
    echo ""
    echo "CITATION REQUIRED:"
    echo "✅ NumPy: Harris et al. (2020) DOI: 10.1038/s41586-020-2649-2"
    
    echo "=========================================="
    
    echo ""
    echo "NEXT STEPS:"
    echo "1. Visualize: ovito ${OUTPUT_XYZ}"
    echo "2. Minimize:  bash scripts/minimize_structure.sh"
    
    if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
        echo "3. Validate grain structure in OVITO:"
        echo "   - Analysis → Grain Segmentation"
        echo "   - Should show ~${GRAINS} grains"
        echo "4. Equilibrate grain boundaries:"
        echo "   - Run NVT at 300K for 10ps"
    elif [[ "${STRUCTURE_TYPE}" == "amorphous" ]]; then
        echo "3. Perform annealing cycle in LAMMPS:"
        echo "   - Heat to 600K, cool to 300K"
        echo "   - This relaxes the amorphous structure"
    fi
    
    echo "=========================================="
else
    echo ""
    echo "❌ ERROR: Output file not created!"
    exit 1
fi