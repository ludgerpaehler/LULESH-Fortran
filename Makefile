FFLAGS = -O3 -flto
CFLAGS = -std=c99 -Wall -O3 -flto
CPPFLAGS =

# Default paths, which can be overwritten by exporting own flags, i.e.
#	`export CLANG=/path/to/own/clang`
# on the command line
LLVM_PATH ?= /home/lpaehler/Work/Dev-Tools/llvm-fortran/f18-llvm-project/build
CLANG ?= $(LLVM_PATH)/build/bin/clang
FLANG ?= $(LLVM_PATH)/build/bin/flang-new

ENZYME_PATH = /home/lpaehler/Work/AutomaticDifferentiation/Enzyme/build/Enzyme/LLVMEnzyme-13.so
LLVM13_PATH = /home/lpaehler/Work/AutomaticDifferentiation/llvm-project/build

default: all
all: lulesh.o

%.fir: %.f90
	# $(FC) $(FFLAGS) $(FFLAGS_add) -c /pwd/$< -o /pwd/$@
	$(LLVM_PATH)/bin/bbc $< -emit-fir --math-runtime=llvm -o $@ -fopenmp

%-unopt.ll: %.fir
	(echo 'target datalayout = "e-m:e-i64:64-f80:128-n8:16:32:64-S128"' && echo 'target triple = "x86_64-pc-linux-gnu"' && $(LLVM_PATH)/bin/tco $^) | $(LLVM_PATH)/bin/opt -S - -o $@
	sed -i -e 's/_QP__enzyme_integer/__enzyme_integer/g' $@ 

combined.ll: lulesh-unopt.ll lulesh_comp_kernels-unopt.ll
	$(LLVM_PATH)/bin/llvm-link lulesh-unopt.ll lulesh_comp_kernels-unopt.ll -S -o $@
	$(LLVM13_PATH)/bin/opt $@ -O2 -o $@ -S

postenzyme.ll: combined.ll
	$(LLVM13_PATH)/bin/opt --enable-new-pm=0 $^ -load=$(ENZYME_PATH) -enzyme -enzyme-loose-types -O2 -o $@ -S
	# $(LLVM13_PATH)/bin/opt $^ -O2 -o $@ -S

lulesh.o: postenzyme.ll
	$(LLVM_PATH)/bin/clang++ -O2 $^ -o $@  -lm -lFortran_main -lFortranRuntime -lFortranDecimal -fopenmp

clean:
	rm -f *.ll *.o

llvm:
	cmake -G Ninja -B llvm-project/build/ -S llvm-project/llvm \
		-DCMAKE_BUILD_TYPE=Release \
		-DLLVM_TARGETS_TO_BUILD="host" \
		-DLLVM_ENABLE_ASSERTIONS=ON \
		-DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra;lld;lldb;flang;openmp" \
		-DLLDB_USE_SYSTEM_DEBUGSERVER=ON \
		-DLLDB_INCLUDE_TESTS=OFF
	cmake --build llvm-project/build/

clean-llvm:
	rm -rf llvm-project/build/
