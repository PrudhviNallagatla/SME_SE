#!/bin/bash
# filepath: /home/rimuru/workspace/scripts/structure_gen.sh
set -euo pipefail

PROJECT_ROOT="/home/rimuru/workspace"
SEED=42

# ============================================
# CONFIGURATION - Edit these variables
# ============================================
# SIZE: Particle diameter in nanometers (nm)
# SHAPE: Available shapes: "sphere", "blob", "ellipsoid"
# GRAINS: 0 = single crystal, >1 = polycrystalline
# AMORPHOUS: 0.0 = crystalline, 0.3 = amorphous
# NI_PERCENT: Nickel percentage (50.0 = equiatomic NiTi)
# ============================================

SIZE=10           # Diameter in nm
SHAPE="blob"      # sphere, blob, or ellipsoid
# either grains or amorphous only - dont use both
GRAINS=0          # 0 = single crystal, >1 = polycrystalline
AMORPHOUS=0.3     # 0.0 = crystalline, 0.3 = amorphous
NI_PERCENT=50.0   # Ni percentage

# Determine structure type for filename
if (( $(python3 -c "print(int(${AMORPHOUS} > 0))") )); then
    TYPE="amorphous"
elif [[ ${GRAINS} -gt 1 ]]; then
    TYPE="polycrystalline"
else
    TYPE="crystalline"
fi

# Output file
OUTPUT="data/structures/niti_${SIZE}nm_${TYPE}.data"

cd "${PROJECT_ROOT}"
mkdir -p data/structures

echo "Generating ${SIZE}nm NiTi nanoparticle (${TYPE})..."

# Build command
CMD="python3 src/niti_np_gen.py --diameter ${SIZE} --shape ${SHAPE} --ni-percent ${NI_PERCENT} --output ${OUTPUT} --seed ${SEED} --lammps"

# Add grains if specified
if [[ ${GRAINS} -gt 1 ]]; then
    CMD="${CMD} --grains ${GRAINS}"
fi

# Add amorphous if specified
if (( $(python3 -c "print(int(${AMORPHOUS} > 0))") )); then
    CMD="${CMD} --amorphous ${AMORPHOUS}"
fi

echo "Command: ${CMD}"
${CMD}

echo "✓ Done: ${OUTPUT}"