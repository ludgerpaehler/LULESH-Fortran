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

To install the Intel oneAPI toolchain, please follow the provided instructions on Intel's [website](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html).
