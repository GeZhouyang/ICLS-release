FC = mpif90

FFLAGS := -O3 -ffixed-line-length-none -fdefault-real-8 
#FFLAGS := -O0 -fbounds-check -g -Wall -Wextra -pedantic -fbacktrace -ffixed-line-length-none -fdefault-real-8 

BIG := # -mcmodel=medium
DBG := -g #-traceback
PROF := # -pg
OMP := #-openmp

# 2decomp&fft
include ../libs/2decomp_fft/src/Makefile.inc
INCLUDE = -I ../libs/2decomp_fft/include

LIB = -L ../libs/2decomp_fft/lib -l2decomp_fft -L ../libs/fft/lib -lfft -llapack

TARGET = ../ghost

SRC = param.f90 common.f90 init.f90 initmpi.f90 initsolver.f90 bound.f90 chkdiv.f90 chkdt.f90 loadd.f90 \
      mom.f90 fillps.f90 zredistribute.f90 solver.f90 correc.f90 output.f90 main.f90 ab2.f90 gfm.f90 \
      levelset.f90 hyperPDE.f90 interface.f90 write_vtk.f90 stopwatch.f90 misc.f90 rib.f90 contact_line.f90

OBJ = $(SRC:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) -o $@ $(OBJ) $(LIB)

.PHONY: clean  
clean:
	rm -f *.o *.mod $(TARGET) .depend *.f90~ *.f90e

veryclean: 
	rm -f *.o *.mod $(TARGET) .depend ../libs/fft/fft.o ../libs/fft/lib/libfft.a;
	cd ../libs/2decomp_fft/src/; make clean; cd ../../../src

libraries:
	cd ../libs/fft; $(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) -c fft.f; ar qc libfft.a fft.o; mv libfft.a lib; cd ../../src;
	cd ../libs/2decomp_fft/src/; make; cd ../../../src

%.o: %.f90
	$(FC) $(INCLUDE) -c -o $@ $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) $<

.depend dep:
	./.makedepo $(SRC) > .depend

include .depend
