#!/bin/bash
set -euo pipefail

# ============================================
# CONFIGURATION - SCIENTIFICALLY VALIDATED
# ============================================
PROJECT_ROOT="/home/rimuru/workspace"
SEED=42

SIZE=20           # Diameter in nm
SHAPE="blob"      # ✅ Use blob for EDM realism
NI_PERCENT=50.0

CORES=-1

# ✅ CRITICAL: Grain size for polycrystalline structures
# Set to 0 for SINGLE CRYSTAL
GRAIN_SIZE=5.0    # Target grain diameter in nm (0 = single crystal)

# ✅ Defect concentrations
VACANCIES=0.01           # 1% bulk vacancies
SURFACE_VACANCIES=0.03   # 3% surface vacancies
ANTISITES=0.005          # 0.5% antisites

# ============================================
# ATOMSK VALIDATION (AUTO-CHECK)
# ============================================
echo "=========================================="
echo "VALIDATING SCIENTIFIC TOOLS"
echo "=========================================="

# Check if Atomsk is installed
if ! command -v atomsk &> /dev/null; then
    echo "❌ CRITICAL ERROR: Atomsk not found!"
    echo ""
    echo "INSTALLATION INSTRUCTIONS:"
    echo "1. Download from: https://atomsk.univ-lille.fr/dl.php"
    echo "2. For Linux:"
    echo "   wget https://atomsk.univ-lille.fr/code/atomsk_b0.13.1_Linux-x86-64.tar.gz"
    echo "   tar -xzf atomsk_b0.13.1_Linux-x86-64.tar.gz"
    echo "   sudo mv atomsk /usr/local/bin/"
    echo "   sudo chmod +x /usr/local/bin/atomsk"
    echo ""
    echo "CITATION REQUIRED:"
    echo "   Hirel, P. (2015). Computer Physics Communications,"
    echo "   197, 212-219. DOI: 10.1016/j.cpc.2015.07.012"
    echo "=========================================="
    exit 1
fi

ATOMSK_VERSION=$(atomsk --version 2>&1 | head -n 1 || echo "unknown")
echo "✅ Atomsk found: ${ATOMSK_VERSION}"

# Validate Python dependencies
if ! python3 -c "import ase, numpy, scipy, psutil" 2>/dev/null; then
    echo "❌ ERROR: Missing Python dependencies"
    echo "   Run: pip install ase numpy scipy tqdm threadpoolctl psutil"
    exit 1
fi
echo "✅ Python dependencies validated"

echo "=========================================="

# ============================================
# STRUCTURE GENERATION
# ============================================

# Determine structure type
if (( $(python3 -c "print(int(${GRAIN_SIZE} > 0))") )); then
    TYPE="nanocrystalline_${GRAIN_SIZE}nm_grains_atomsk"
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
    echo "Grain size:    ${GRAIN_SIZE} nm (Atomsk polycrystal)"
    echo "Method:        Atomsk (MANDATORY - peer-reviewed)"
else
    echo "Method:        Single crystal (no grain boundaries)"
fi

echo "Seed:          ${SEED}"
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

# ✅ CRITICAL: Add grain structure ONLY if polycrystalline
if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
    CMD="${CMD} --grain-size ${GRAIN_SIZE}"
    CMD="${CMD} --use-atomsk"  # ✅ MANDATORY for polycrystals
fi

# Add defects
CMD="${CMD} --vacancies ${VACANCIES}"
CMD="${CMD} --surface-vacancies ${SURFACE_VACANCIES}"
CMD="${CMD} --antisites ${ANTISITES}"

echo "Executing: ${CMD}"
echo "=========================================="

# Execute
eval ${CMD}

# ✅ VALIDATION: Check that Atomsk was actually used
if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
    if [[ -f "${OUTPUT}" ]]; then
        # Check if file contains Atomsk signature
        if ! grep -q "Atomsk" "${OUTPUT}" 2>/dev/null; then
            echo ""
            echo "⚠️  WARNING: Output file doesn't contain Atomsk signature!"
            echo "   Atomsk may have failed silently - check logs above"
            echo ""
        else
            echo ""
            echo "✅ ATOMSK VALIDATION PASSED"
        fi
    fi
fi

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
    
    if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
        echo ""
        echo "CITATION REQUIRED:"
        echo "✅ Atomsk: Hirel (2015) DOI: 10.1016/j.cpc.2015.07.012"
        echo "✅ ASE: Larsen et al. (2017) DOI: 10.1088/1361-648X/aa680e"
    else
        echo ""
        echo "CITATION REQUIRED:"
        echo "✅ ASE: Larsen et al. (2017) DOI: 10.1088/1361-648X/aa680e"
    fi
    
    echo "=========================================="
    
    echo ""
    echo "NEXT STEPS:"
    echo "1. Visualize: ovito ${OUTPUT_XYZ}"
    echo "2. Minimize:  bash scripts/minimize_structure.sh"
    
    if [[ "${STRUCTURE_TYPE}" == "polycrystalline" ]]; then
        echo "3. Validate grain boundaries in OVITO:"
        echo "   - Analysis → Grain Segmentation"
        echo "   - Should show ~$(python3 -c "print(int((${SIZE}/${GRAIN_SIZE})**3))") grains"
        echo "   - GB thickness should be 2-4 Å"
    fi
    
    echo "=========================================="
else
    echo ""
    echo "❌ ERROR: Output file not created!"
    exit 1
fi