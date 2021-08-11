FFLAGS = -O3 -flto
CFLAGS = -std=c99 -Wall -O3 -flto
CPPFLAGS =

LOAD = /path

CLANG = /home/lpaehler/Work/Dev-Tools/llvm-fortran/f18-llvm-project/build/bin/clang

LLVM_PATH = /home/lpaehler/Work/Dev-Tools/llvm-fortran/f18-llvm-project/build

default: all
all: lulesh.o

%-unopt.ll: %.f90
	# $(FC) $(FFLAGS) $(FFLAGS_add) -c /pwd/$< -o /pwd/$@
	(echo 'target datalayout = "e-m:e-i64:64-f80:128-n8:16:32:64-S128"' && echo 'target triple = "x86_64-pc-linux-gnu"' && $(LLVM_PATH)/bin/bbc $< -emit-fir -o - --math-runtime=llvm | $(LLVM_PATH)/bin/tco) | $(LLVM_PATH)/bin/opt - -o $@

%-raw.ll: %-unopt.ll
	# opt $^ -load=$(LOAD) -enzyme -o $@ -S
	opt $^ -o $@ -S

%-opt.ll: %-raw.ll
	opt $^ -o $@ -S

lulesh.o: lulesh-opt.ll
	$(LLVM_PATH)/bin/clang++ -O2 $^ -o $@  -lm -lFortran_main -lFortranRuntime -lFortranDecimal

# FLAGS: -cpp -fallow-invalid-boz -ffree-line-length-none 