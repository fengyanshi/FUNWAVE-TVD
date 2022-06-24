CPP = /usr/bin/cpp
CPPFLAGS = -DINTEL -P -traditional
INCS        = $(IOINCS) $(GOTMINCS)
LIBS        = $(PV3LIB) $(CLIB)  $(PARLIB) $(IOLIBS) $(MPILIB) $(GOTMLIB)

ifdef DEBUG
  DEBFLGS  = -g -check -warn 
  OPT      = -O0
else
  DEBFLGS  =
  OPT      = -O2
endif

ifdef SEQU
  FC=/opt/intel/oneapi/compiler/2022.1.0/linux/bin/intel64/ifort
  FFLAGS=-fPIC
  LFLAGS=-shared-intel
else
  FC = /opt/intel/oneapi/mpi/2021.6.0/bin/mpiifort
  FFLAGS=-fPIC
  LFLAGS=-shared-intel
  #/opt/intel/oneapi/mpi/2021.6.0/bin/mpirun -np 8 ./funwave
endif

ifdef DOUBLE_PRECISION
  FFLAGS += -fr8
endif

