#!/bin/bash
# filepath: d:\MD Sims\SME_SE\all_dependencies.sh

###############################################################################
# LAMMPS Full Installation Script for WSL/Docker/DevContainers
# Optimized for CLI-only environments (no GUI)
# 
# Usage: 
#   chmod +x all_dependencies.sh
#   ./all_dependencies.sh
#
# For Docker:
#   docker run -it -v $(pwd):/workspace ubuntu:22.04
#   cd /workspace && bash all_dependencies.sh
#
# Author: PrudhviNallagatla
###############################################################################

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Color output for better visibility
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log() { echo -e "${GREEN}[$(date +'%H:%M:%S')]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

# Detect environment
if grep -qi microsoft /proc/version 2>/dev/null; then
    ENV_TYPE="WSL"
elif [ -f /.dockerenv ]; then
    ENV_TYPE="Docker"
else
    ENV_TYPE="Linux"
fi

log "Detected environment: ${ENV_TYPE}"

# Detect number of cores for parallel compilation
NPROC=$(nproc)
log "Detected ${NPROC} CPU cores"

# Set non-interactive for apt (critical for Docker)
export DEBIAN_FRONTEND=noninteractive
export TZ=Etc/UTC

#==============================================================================
# STEP 1: Remove old LAMMPS installation
#==============================================================================
log "[1/8] Checking for existing LAMMPS installation..."

if command -v lmp &> /dev/null || command -v lammps &> /dev/null; then
    warn "Found existing LAMMPS installation. Removing..."
    
    sudo rm -f /usr/local/bin/lmp /usr/local/bin/lammps 2>/dev/null || true
    sudo rm -rf /usr/local/share/lammps 2>/dev/null || true
    sudo rm -rf /usr/local/lib/liblammps* 2>/dev/null || true
    sudo rm -rf /tmp/lammps 2>/dev/null || true
    
    # Remove Python LAMMPS interface
    pip3 uninstall -y lammps 2>/dev/null || true
    
    log "✓ Old LAMMPS removed"
else
    log "No existing LAMMPS found (clean install)"
fi

#==============================================================================
# STEP 2: Update system
#==============================================================================
log "[2/8] Updating system packages..."
sudo apt-get update -qq
sudo apt-get upgrade -yqq

#==============================================================================
# STEP 3: Install system dependencies (NO GUI tools)
#==============================================================================
log "[3/8] Installing build dependencies and CLI tools..."
sudo apt-get install -yqq --no-install-recommends \
    build-essential \
    cmake \
    git \
    git-lfs \
    ca-certificates \
    wget \
    curl \
    vim \
    nano \
    tmux \
    tree \
    unzip \
    zip \
    pigz \
    zstd \
    xz-utils \
    rsync \
    jq \
    gfortran \
    g++ \
    gcc \
    make \
    pkg-config \
    libopenmpi-dev \
    openmpi-bin \
    libblas-dev \
    liblapack-dev \
    libfftw3-dev \
    libeigen3-dev \
    libpng-dev \
    libjpeg-dev \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    libhdf5-dev \
    libnetcdf-dev \
    zlib1g-dev \
    libgsl-dev \
    graphviz \
    imagemagick \
    gnuplot-nox \
    ffmpeg \
    bc \
    time

log "✓ System dependencies installed"

#==============================================================================
# STEP 4: Install Python packages
#==============================================================================
log "[4/8] Installing Python scientific packages..."

# Upgrade pip first
pip3 install --break-system-packages --no-cache-dir --upgrade pip setuptools wheel

# Install packages in batches to avoid memory issues
log "  Installing core scientific stack..."
pip3 install --break-system-packages --no-cache-dir \
    numpy \
    scipy \
    matplotlib \
    pandas

log "  Installing molecular dynamics tools..."
pip3 install --break-system-packages --no-cache-dir \
    ase \
    ovito \
    MDAnalysis \
    pymatgen \
    freud-analysis

log "  Installing analysis and utilities..."
pip3 install --break-system-packages --no-cache-dir \
    jupyter \
    jupyterlab \
    ipykernel \
    tqdm \
    pyyaml \
    h5py \
    netcdf4 \
    seaborn \
    plotly \
    networkx \
    scikit-learn

log "  Installing advanced tools..."
pip3 install --break-system-packages --no-cache-dir \
    lmfit \
    statsmodels \
    sympy \
    numba \
    cython \
    atomman \
    matscipy \
    phonopy \
    spglib \
    mendeleev \
    periodictable \
    tabulate \
    rich \
    click \
    loguru

log "✓ Python packages installed"

#==============================================================================
# STEP 5: Install Atomsk (CRITICAL for polycrystalline structures)
#==============================================================================
# ⚠️  ATOMSK MUST BE INSTALLED MANUALLY ⚠️
#
# Atomsk is ESSENTIAL for creating grain boundaries and polycrystalline structures
# The automated download often fails due to server issues.
#
# MANUAL INSTALLATION STEPS:
# --------------------------
# 1. Visit: https://atomsk.univ-lille.fr/install.php
# 
# 2. Download the latest Linux binary:
#    wget https://atomsk.univ-lille.fr/code/atomsk_b0.13.1_Linux-x86-64.tar.gz
#    (or use your browser if wget fails)
#
# 3. Extract and install:
#    tar -xzf atomsk_b0.13.1_Linux-x86-64.tar.gz
#    sudo mv atomsk /usr/local/bin/
#    sudo chmod +x /usr/local/bin/atomsk
#
# 4. Verify installation:
#    atomsk --version
#    atomsk --help
#
# 5. Test polycrystal creation:
#    atomsk --create fcc 4.05 Al aluminum.xsf
#    atomsk --polycrystal aluminum.xsf polycrystal.cfg -wrap
#
# ALTERNATIVE: Install via Conda (if available):
#    conda install -c conda-forge atomsk
#
# ALTERNATIVE: Compile from source:
#    git clone https://github.com/pierrehirel/atomsk.git
#    cd atomsk/src
#    make atomsk
#    sudo cp atomsk /usr/local/bin/
#
# For Docker: Add this to your Dockerfile AFTER this script runs:
#    RUN wget https://atomsk.univ-lille.fr/code/atomsk_b0.13.1_Linux-x86-64.tar.gz && \
#        tar -xzf atomsk_b0.13.1_Linux-x86-64.tar.gz && \
#        mv atomsk /usr/local/bin/ && \
#        chmod +x /usr/local/bin/atomsk && \
#        rm atomsk_b0.13.1_Linux-x86-64.tar.gz
#
#==============================================================================

log "[5/8] Atomsk installation SKIPPED - requires manual installation"
warn "═══════════════════════════════════════════════════════════════"
warn "  ⚠️  ATOMSK MUST BE INSTALLED MANUALLY FOR POLYCRYSTALS  ⚠️"
warn "═══════════════════════════════════════════════════════════════"
warn ""
warn "Quick install (run after this script completes):"
warn ""
warn "  cd /tmp"
warn "  wget https://atomsk.univ-lille.fr/code/atomsk_b0.13.1_Linux-x86-64.tar.gz"
warn "  tar -xzf atomsk_b0.13.1_Linux-x86-64.tar.gz"
warn "  sudo mv atomsk /usr/local/bin/"
warn "  sudo chmod +x /usr/local/bin/atomsk"
warn "  atomsk --version"
warn ""
warn "Official docs: https://atomsk.univ-lille.fr/doc/en/"
warn "═══════════════════════════════════════════════════════════════"
warn ""

# Check if Atomsk already exists (from previous manual install)
if command -v atomsk &> /dev/null; then
    log "✓ Atomsk already installed: $(atomsk --version 2>&1 | head -1)"
else
    warn "Atomsk NOT found - remember to install manually!"
fi

#==============================================================================
# STEP 6: Download LAMMPS
#==============================================================================
log "[6/8] Downloading LAMMPS from GitHub..."
cd /tmp

# Remove any existing LAMMPS source
rm -rf lammps

# Clone stable version (faster than full history)
if ! git clone --depth 1 --branch stable_2Aug2023_update3 https://github.com/lammps/lammps.git; then
    error "Failed to download LAMMPS from GitHub"
fi

cd lammps
log "✓ LAMMPS source downloaded ($(git describe --tags))"

#==============================================================================
# STEP 7: Compile LAMMPS with ALL packages
#==============================================================================
log "[7/8] Compiling LAMMPS with ALL packages..."
log "  This will take 15-30 minutes depending on CPU..."

mkdir -p build
cd build

# Configure with CMake
cmake ../cmake \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_MPI=yes \
    -D BUILD_OMP=yes \
    -D BUILD_SHARED_LIBS=yes \
    -D CMAKE_INSTALL_PREFIX=/usr/local \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_C_COMPILER=mpicc \
    -D CMAKE_Fortran_COMPILER=mpifort \
    -D PKG_ASPHERE=yes \
    -D PKG_BODY=yes \
    -D PKG_CLASS2=yes \
    -D PKG_COLLOID=yes \
    -D PKG_COMPRESS=yes \
    -D PKG_CORESHELL=yes \
    -D PKG_DIPOLE=yes \
    -D PKG_GRANULAR=yes \
    -D PKG_KSPACE=yes \
    -D PKG_MANYBODY=yes \
    -D PKG_MC=yes \
    -D PKG_MEAM=yes \
    -D PKG_MISC=yes \
    -D PKG_MOLECULE=yes \
    -D PKG_OPT=yes \
    -D PKG_PERI=yes \
    -D PKG_QEQ=yes \
    -D PKG_REPLICA=yes \
    -D PKG_RIGID=yes \
    -D PKG_SHOCK=yes \
    -D PKG_SNAP=yes \
    -D PKG_SRD=yes \
    -D PKG_REAXFF=yes \
    -D PKG_SMTBQ=yes \
    -D PKG_PYTHON=yes \
    -D PKG_EXTRA-COMPUTE=yes \
    -D PKG_EXTRA-DUMP=yes \
    -D PKG_EXTRA-FIX=yes \
    -D PKG_EXTRA-MOLECULE=yes \
    -D PKG_EXTRA-PAIR=yes \
    -D PKG_ML-SNAP=yes \
    -D PKG_PHONON=yes \
    -D PKG_ATC=yes \
    -D PKG_AWPMD=yes \
    -D PKG_COLVARS=yes \
    -D PKG_DIFFRACTION=yes \
    -D PKG_EFF=yes \
    -D PKG_FEP=yes \
    -D PKG_H5MD=yes \
    -D PKG_NETCDF=yes \
    -D PKG_OPENMP=yes \
    -D PKG_PTM=yes \
    -D PKG_QTB=yes \
    -D PKG_TALLY=yes \
    -D PKG_VORONOI=yes \
    -D PKG_USER-CGDNA=yes \
    -D PKG_USER-DRUDE=yes \
    -D PKG_USER-MESONT=yes \
    -D PKG_USER-PTM=yes \
    > /tmp/lammps_cmake.log 2>&1

if [ $? -ne 0 ]; then
    error "CMake configuration failed. Check /tmp/lammps_cmake.log"
fi

log "  CMake configuration complete"
log "  Building with ${NPROC} parallel jobs..."

# Build (with progress indicator)
if ! make -j${NPROC} > /tmp/lammps_build.log 2>&1; then
    error "LAMMPS compilation failed. Check /tmp/lammps_build.log"
fi

log "✓ LAMMPS compiled successfully"

# Install
log "  Installing LAMMPS to /usr/local..."
if ! sudo make install > /tmp/lammps_install.log 2>&1; then
    error "LAMMPS installation failed. Check /tmp/lammps_install.log"
fi

log "✓ LAMMPS installed"

#==============================================================================
# STEP 8: Install LAMMPS Python interface
#==============================================================================
log "[8/8] Installing LAMMPS Python interface..."
cd ../python

if ! pip3 install --break-system-packages --no-cache-dir . > /tmp/lammps_python.log 2>&1; then
    warn "LAMMPS Python interface installation failed (check /tmp/lammps_python.log)"
else
    log "✓ LAMMPS Python interface installed"
fi

#==============================================================================
# VERIFICATION
#==============================================================================
log "Verifying installation..."

# Check LAMMPS binary
if ! command -v lmp &> /dev/null; then
    error "LAMMPS binary not found in PATH"
fi

# Check version
LAMMPS_VERSION=$(lmp -help 2>&1 | grep "LAMMPS" | head -1)
log "✓ ${LAMMPS_VERSION}"

# Check Python interface
if python3 -c "from lammps import lammps; lammps()" &> /dev/null; then
    log "✓ LAMMPS Python interface working"
else
    warn "LAMMPS Python interface not working (non-critical)"
fi

# Check key packages
log "Checking enabled packages..."
lmp -help 2>&1 | grep -A50 "Installed packages" | grep -E "MEAM|VORONOI|PYTHON|EXTRA-DUMP" || warn "Some packages missing"

#==============================================================================
# CLEANUP
#==============================================================================
log "Cleaning up temporary files..."
cd ~
sudo rm -rf /tmp/lammps
sudo apt-get autoremove -yqq
sudo apt-get clean

#==============================================================================
# POST-INSTALL CONFIGURATION
#==============================================================================
log "Setting up environment..."

# Create .bashrc additions (only if not exists)
if ! grep -q "# LAMMPS aliases" ~/.bashrc 2>/dev/null; then
    cat >> ~/.bashrc << 'EOF'

# LAMMPS aliases
alias lammps='lmp'
alias lmp-serial='lmp'
alias lmp-mpi='mpirun -np $(nproc) lmp'
alias lmp-omp='lmp -sf omp -pk omp $(nproc)'

# Python aliases
alias py='python3'
alias ipy='ipython'

# Analysis aliases
alias ovito='ovito'
EOF
    log "✓ Aliases added to ~/.bashrc"
else
    log "✓ Aliases already exist in ~/.bashrc"
fi

#==============================================================================
# FINAL REPORT
#==============================================================================
echo ""
echo "=========================================="
echo -e "${GREEN}✓ INSTALLATION COMPLETE${NC}"
echo "=========================================="
echo "Environment: ${ENV_TYPE}"
echo "LAMMPS:      $(lmp -help 2>&1 | grep "LAMMPS" | head -1)"
echo "Python:      $(python3 --version)"
echo "Cores:       ${NPROC}"
echo "=========================================="
echo ""

# Check Atomsk status
if command -v atomsk &> /dev/null; then
    echo -e "${GREEN}✓ Atomsk:     $(atomsk --version 2>&1 | head -1)${NC}"
else
    echo -e "${RED}✗ Atomsk:     NOT INSTALLED${NC}"
    echo ""
    echo -e "${YELLOW}⚠️  INSTALL ATOMSK MANUALLY FOR POLYCRYSTAL WORK:${NC}"
    echo ""
    echo "  cd /tmp"
    echo "  wget https://atomsk.univ-lille.fr/code/atomsk_b0.13.1_Linux-x86-64.tar.gz"
    echo "  tar -xzf atomsk_b0.13.1_Linux-x86-64.tar.gz"
    echo "  sudo mv atomsk /usr/local/bin/"
    echo "  sudo chmod +x /usr/local/bin/atomsk"
    echo ""
fi

echo "=========================================="
echo ""
echo "Quick start:"
echo "  lmp -in your_script.lmp"
echo "  mpirun -np ${NPROC} lmp -in your_script.lmp"
echo "  python3 your_analysis.py"
echo ""
echo "Installed Python packages:"
pip3 list 2>/dev/null | grep -E "numpy|scipy|matplotlib|ovito|ase|MDAnalysis|freud" || echo "  (list error, but packages installed)"
echo ""
echo "LAMMPS packages:"
lmp -help 2>&1 | grep -A50 "Installed packages" | head -20
echo ""
if [ "${ENV_TYPE}" = "WSL" ]; then
    echo "WSL Note: Restart your shell or run: source ~/.bashrc"
fi
echo "=========================================="