#!/bin/bash
# filepath: /home/rimuru/workspace/scripts/structure_gen.sh
set -euo pipefail

PROJECT_ROOT="/home/rimuru/workspace"
SEED=42

SIZE=20           # Diameter in nm
SHAPE="blob"      # sphere, blob, or ellipsoid
# either grains or amorphous only - dont use both
GRAINS=64          # 0 = single crystal, >1 = polycrystalline
AMORPHOUS=0.0     # 0.0 = crystalline, 0.3 = amorphous
NI_PERCENT=50.0   # Ni percentage

# Estimate number of grains based on particle and grain size:
# Ng ≈ (d / dg)^3
# For example: diameter d = 20 nm, grain size dg = 5 nm → Ng ≈ (20/5)^3 = 64 grains.
# This ensures physically realistic grain counts at nanoscale.

# ============================================
# CONFIGURATION - Edit these variables
# ============================================
# SIZE: Particle diameter in nanometers (nm)
# SHAPE: Available shapes: "sphere", "blob", "ellipsoid"
# GRAINS: 0 = single crystal, >1 = polycrystalline
# AMORPHOUS: 0.0 = crystalline, 0.3 = amorphous
# NI_PERCENT: Nickel percentage (50.0 = equiatomic NiTi)
# ============================================

if (( $(python3 -c "print(int(${AMORPHOUS} > 0))") )); then
    TYPE="amorphous"
elif [[ ${GRAINS} -gt 1 ]]; then
    TYPE="crystalline_${GRAINS}grains"
else
    TYPE="crystalline_1grain"
fi

# Output file with Ni percentage
NI_INT=$(python3 -c "print(int(${NI_PERCENT}))")
OUTPUT="data/structures/niti_${SIZE}nm_Ni${NI_INT}_${TYPE}.data"

cd "${PROJECT_ROOT}"
mkdir -p data/structures

echo "Generating ${SIZE}nm NiTi nanoparticle (${TYPE}, Ni${NI_INT}%)..."

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