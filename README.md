# LULESH Fortran

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

