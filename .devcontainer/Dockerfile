# Bare-bones Ubuntu container
# All installations happen via install_lammps.sh

FROM ubuntu:22.04

LABEL maintainer="PrudhviNallagatla"
LABEL description="Minimal Ubuntu base for LAMMPS installation"

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Kolkata

# Create user
ARG USERNAME=rimuru
ARG USER_UID=1000
ARG USER_GID=1000

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    && apt-get update \
    && apt-get install -y sudo \
    && echo "$USERNAME ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME \
    && rm -rf /var/lib/apt/lists/*

# Create workspace
RUN mkdir -p /home/$USERNAME/workspace && chown -R $USERNAME:$USERNAME /home/$USERNAME

USER $USERNAME
WORKDIR /home/$USERNAME/workspace

CMD ["/bin/bash"]