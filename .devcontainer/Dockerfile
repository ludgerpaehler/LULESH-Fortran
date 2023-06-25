FROM ubuntu:22.04

RUN apt-get update && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    curl ca-certificates gpg && \
  rm -rf /var/lib/apt/lists/*

# repository to install Intel(R) oneAPI Libraries
RUN curl -fsSL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB | gpg --dearmor - | tee /usr/share/keyrings/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.gpg >/dev/null
RUN echo "deb [signed-by=/usr/share/keyrings/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.gpg] https://apt.repos.intel.com/oneapi all main" > /etc/apt/sources.list.d/oneAPI.list

RUN apt-get update && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    curl ca-certificates gpg-agent software-properties-common openssh-client \
    build-essential pkg-config \
    git ninja-build make \
    intel-oneapi-compiler-fortran gfortran \

ENV LANG=C.UTF-8