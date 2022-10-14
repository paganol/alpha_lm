
HOST=$(shell hostname)

ifneq ($(findstring cori,$(HOST)), )
#Cori
F90C=ftn
F90FLAGS=-O3 -march=native -Wall -fopenmp -x f95-cpp-input -ffree-line-length-0 #-fcheck=bounds#-ffree-form 
INCLUDES=-I$(cfitsio_PREFIX)/include -I$(healpix_PREFIX)/include $(shell toastconfig --cppflags) -I$(fftw_PREFIX)/include
LIBS=-L$(healpix_PREFIX)/lib -L$(cfitsio_PREFIX)/lib -lhealpix -lchealpix -lcfitsio -lm -lgfortran -L$(getdata_PREFIX)/lib -lgetdata -L$(fftw_PREFIX)/lib -lfftw3
endif

ifneq ($(findstring log0,$(HOST)), )
#Magique4
F90C= mpif90
F90FLAGS= -mkl -qopenmp -O3 -no-prec-div -heap-arrays 524288 #2097152 -fpp #-no-heap-arrays #-heap-arrays #30000 #2097152 #524288 #30000 #16777216 #2097152 #16777216 #2097152 #524288#-fpp -check all -check bounds 
INCLUDES= -I. -I$(HEALPIX)/include -I$(MKL_MOD)
LIBS= -L. -L$(HEALPIX)/lib -lhealpix -L/wrk/lpagano/likelihoods/plc-2.0/lib -lcfitsio -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
endif


ifneq ($(findstring coka,$(HOST)), )
F90C= mpiifort
F90FLAGS= -O3 -march=native -funroll-loops -fopenmp -check bounds #-g -check all -fpe0 -warn -traceback -debug extended -check bounds
INCLUDES= -I. -I/usr/local/include -I$(HEALPIX)/include/ -I/usr/include
LIBS= -L. -L$(HEALPIX)/lib -lhealpix -lcfitsio -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
endif

ifneq ($(findstring r000u,$(HOST)), )
#marconi - gcc
#F90C    = mpif90
###F90FLAGS = -O3 -fopenmp -ffree-form -x f95-cpp-input -ffast-math -march=native -funroll-loops -fno-second-underscore -fPIC 
#F90FLAGS = -O3 -fPIC -ffree-line-length-none -fopenmp -Wall -fcheck=bounds -fbacktrace
#INCLUDES = -I. -I$(HEALPIX)/include -I$(CFITSIO)/include
#LIBS = -L$(HEALPIX)/lib -lhealpix -lsharp -L$(CFITSIO)/lib -lcfitsio
#marconi - intel
F90C= mpiifort
#https://www.spec.org/hpc2021/flags/Intel_compiler_flags.2021-10-20.html
#https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compiler-reference/compiler-options/code-generation-options/march.htm
F90FLAGS= -Ofast -mkl -march=skylake -funroll-loops -qopenmp -no-prec-div #-check all -traceback -warn#-no-prec-div -check bounds -traceback #-g -check all -fpe0 -warn -traceback -debug extended -check bounds
INCLUDES= -I. -I/usr/local/include -I$(HEALPIX_INCLUDE) -I/usr/include
LIBS= -L. -L$(HEALPIX_LIB) -lhealpix -L$(CFITSIO_LIB) -lcfitsio -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
endif

ifneq ($(findstring local,$(HOST)), )
#my macbook
F90C    = mpif90 
#F90FLAGS = -O3 -fopenmp -ffree-form -x f95-cpp-input -ffast-math -march=native -funroll-loops -fno-second-underscore -fPIC 
F90FLAGS = -O3 -fPIC -ffree-line-length-none -fopenmp -Wall -fcheck=bounds -fbacktrace -x f95-cpp-input -march=native
INCLUDES = -I. -I/usr/local/include -I$(HEALPIX)/include -I/usr/local/Cellar/cfitsio/4.1.0/include -I/usr/local/Cellar/shtools/4.10/include
LIBS = -L$(HEALPIX)/lib -lhealpix  -L/usr/local/Cellar/cfitsio/4.1.0/lib -lcfitsio
endif

OBJFILES = d1mach.o settings.o driver.o utils.o wigner.o wigner1.o EB_estimator.o 

all : EB_estimator


%.o: %.f90
	$(F90C) $(F90FLAGS) $(INCLUDES) -c $*.f90

%.o: %.f
	$(F90C) $(F90FLAGS) $(INCLUDES) -c $*.f

EB_estimator : $(OBJFILES)
	$(F90C) $(OBJFILES) -o ../EB_estimator $(LIBS)

clean :
	rm -f *.o *.mod ../EB_estimator
