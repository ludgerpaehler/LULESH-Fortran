ARG UBUNTU_VERSION=22.04

FROM ubuntu:$UBUNTU_VERSION as base

#ARG USERNAME=vscode
#ARG USER_UID=1000
#ARG USER_GID=$USER_UID

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates

RUN apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common curl gpg gnupg2 openssh-client build-essential pkg-config git make git ssh gcc g++ gfortran wget

# Download the oneAPI key to the system keyring
RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor > /usr/share/keyrings/oneapi-archive-keyring.gpg

# Add signed entry to apt sources and configure the APT client to use the Intel repository
RUN echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" > /etc/apt/sources.list.d/oneAPI.list

# Install oneAPI itself
RUN apt-get update \
    && apt-get install -y intel-oneapi-compiler-dpcpp-cpp-and-cpp-classic intel-oneapi-compiler-fortran

# Cleanup
RUN apt-get autoremove -y --purge \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*