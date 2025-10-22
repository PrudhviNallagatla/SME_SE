#!/bin/bash

###############################################################################
# LAMMPS Full Installation Script for Ubuntu
# For Azure VM: 128GB RAM, 8 cores
# Also works on laptop with 8GB RAM (just takes longer to compile)
#
# Usage: 
#   chmod +x install_lammps.sh
#   ./install_lammps.sh
#
# Author: PrudhviNallagatla
###############################################################################

set -e  # Exit on any error

# Detect number of cores for parallel compilation
NPROC=$(nproc)
echo "Detected ${NPROC} CPU cores"

# Update system
echo ""
echo "[1/6] Updating system packages..."
sudo apt-get update
sudo apt-get upgrade -y

# Install all dependencies
echo ""
echo "[2/6] Installing build dependencies and tools..."
sudo apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    ca-certificates \
    wget \
    curl \
    vim \
    nano \
    htop \
    tmux \
    tree \
    unzip \
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
    ffmpeg \
    python3 \
    python3-dev \
    python3-pip \
    python3-numpy \
    python3-scipy \
    python3-matplotlib \
    libhdf5-dev \
    libnetcdf-dev \
    zlib1g-dev \
    libgsl-dev \
    libvtk9-dev

# Install Python packages
echo ""
echo "[3/6] Installing Python scientific packages..."
pip3 install --no-cache-dir --upgrade pip
pip3 install --no-cache-dir \
    numpy \
    scipy \
    matplotlib \
    pandas \
    ovito==3.9.1 \
    ase \
    jupyter \
    jupyterlab \
    ipykernel \
    tqdm \
    pyyaml \
    h5py \
    seaborn \
    scikit-learn \
    MDAnalysis \
    pymatgen

# Download LAMMPS
echo ""
echo "[4/6] Downloading LAMMPS from GitHub..."
cd /tmp
if [ -d "lammps" ]; then
    rm -rf lammps
fi
git clone --depth 1 --branch stable_2Aug2023_update3 https://github.com/lammps/lammps.git
cd lammps

# Build LAMMPS with ALL packages
echo ""
echo "[5/6] Compiling LAMMPS with ALL packages (this will take 15-30 minutes)..."
mkdir -p build
cd build

cmake ../cmake \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_MPI=yes \
    -D BUILD_OMP=yes \
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
    -D PKG_VTK=yes

echo ""
echo "Building with ${NPROC} parallel jobs..."
make -j${NPROC}

echo ""
echo "Installing LAMMPS to /usr/local..."
sudo make install

# Verify installation
echo ""
echo "[6/6] Verifying installation..."
if command -v lmp &> /dev/null; then
    echo "✓ LAMMPS installed successfully!"
    lmp -help | head -20
else
    echo "✗ LAMMPS installation failed!"
    exit 1
fi

# Cleanup
echo ""
echo "Cleaning up temporary files..."
cd ~
sudo rm -rf /tmp/lammps
sudo apt-get autoremove -y
sudo apt-get clean

echo ""
echo "=================================================="
echo "Installation completed successfully!"
echo "Finished at: $(date)"
echo "=================================================="
echo ""
echo "Quick start:"
echo "  - Run LAMMPS: lmp -in your_script.in"
echo "  - Run parallel: mpirun -np 4 lmp -in your_script.in"
echo "  - Python: python3"
echo "  - Jupyter: jupyter lab --ip=0.0.0.0 --no-browser"
echo ""
echo "Python packages installed:"
pip3 list | grep -E "numpy|scipy|matplotlib|ovito|ase|pandas|jupyter"
echo ""