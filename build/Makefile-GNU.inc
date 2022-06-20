CPP = /usr/bin/cpp
CPPFLAGS = -P -traditional
INCS        = $(IOINCS) $(GOTMINCS)
LIBS        = $(PV3LIB) $(CLIB)  $(PARLIB) $(IOLIBS) $(MPILIB) $(GOTMLIB)

ifdef DEBUG
  # -Wextra -ffpe-trap=underflow,denormal
  DEBFLGS  = -g -fcheck=all -Wall -Werror -fmax-errors=1 -ffpe-trap=invalid,zero,overflow -O0
else
  DEBFLGS  = -O3
endif

ifdef SEQU
  FC=gfortran
  FFLAGS=-nocpp -fmax-stack-var-size=80000 -fimplicit-none -fPIC
  LFLAGS=
else
  FC = mpif90
  FFLAGS=-nocpp -fmax-stack-var-size=80000 -fimplicit-none -fPIC
  LFLAGS=
endif

ifdef FLAG_PRECISION
  FFLAGS += -fmax-stack-var-size=80000 -fdefault-real-8
else
  FFLAGS += -fmax-stack-var-size=40000
endif

