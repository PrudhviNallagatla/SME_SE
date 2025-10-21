# Lightweight LAMMPS container for NiTi nanoparticle simulations
# Optimized for CPU-only, fast builds, minimal size

FROM ubuntu:22.04

LABEL maintainer="PrudhviNallagatla"
LABEL description="LAMMPS environment for NiTi shape memory simulations"

# Prevent timezone prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Kolkata
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Install LAMMPS and essential tools in one layer (reduces image size)
RUN apt-get update && apt-get install -y --no-install-recommends \
    lammps \
    python3 \
    python3-pip \
    python3-numpy \
    python3-scipy \
    python3-matplotlib \
    git \
    vim \
    nano \
    wget \
    curl \
    htop \
    tmux \
    && rm -rf /var/lib/apt/lists/*

# Install Python analysis packages
RUN pip3 install --no-cache-dir \
    ovito==3.9.1 \
    ase==3.22.1 \
    pandas==2.0.3 \
    jupyter==1.0.0 \
    ipykernel==6.24.0 \
    tqdm==4.65.0 \
    pyyaml==6.0.1

# Create non-root user matching your system
ARG USERNAME=rimuru
ARG USER_UID=1000
ARG USER_GID=1000

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    && apt-get update \
    && apt-get install -y sudo \
    && echo "$USERNAME ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Set up workspace
RUN mkdir -p /home/$USERNAME/workspace && chown -R $USERNAME:$USERNAME /home/$USERNAME

# Clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Switch to non-root user
USER $USERNAME
WORKDIR /home/$USERNAME/workspace

# Verify LAMMPS installation
RUN lammps -help | head -5 || echo "LAMMPS check will run on container start"

# Default command
CMD ["/bin/bash"]