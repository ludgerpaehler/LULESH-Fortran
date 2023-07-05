# LULESH Fortran

## Compilation

tbd

## Building the Compiler Toolchain

### LLVM Toolchain

For the llvm-project dependency, we can either clone recursively or can initialize the submodules post-hoc

```bash
git clone --recursive https://github.com/ludgerpaehler/LULESH-Fortran.git
```

or for initializing the submodule

```bash
git submodule update --init --recursive
```

To then build the LLVM toolchain one just has to call onto the Makefile from the root of the repository

```bash
make llvm
```

### Intel OneAPI Toolchain

To install the Intel oneAPI toolchain, please follow the provided instructions on Intel's [website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html). To be able to compile LULESH with Intel compilers you require the following two components from Intel oneAPI, which have to be installed in the specified order:

1. [Intel oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=mac&distributions=online)
2. [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?operatingsystem=mac&distributions=online)
