# LULESH Fortran

## LLVM Production Build

Used FORTRAN compiler

```bash
git clone https://github.com/wsmoses/f18-llvm-project/tree/fir-dev
```

Go into the repo, create a build directory, and 

```bash
cd f18-llvm-project
mkdir build
cd build
```

And then build F18

```bash
cmake -G Ninja ../llvm/ -DLLVM_TARGETS_TO_BUILD="host" -DLLVM_ENABLE_PROJECTS="clang;flang;parallel-libs;openmp" -DLLVM_ENABLE_PLUGINS=ON -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_ASSERTIONS=ON
ninja
```


## GFortran Debug Build

Using GNU Fortran (GCC) 10.2.1 packaged with Fedora

```bash
gfortran lulesh.f90 lulesh_comp_kernels.f90 -fallow-invalid-boz -cpp
```

