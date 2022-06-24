SRC=../src
DST=../build
# FLAG_PRECISION = -DDOUBLE_PRECISION 
FLAG_CARTESIAN = -DCARTESIAN
# FLAG_1  = -DCOUPLING
# FLAG_2  = -DZALPHA
# FLAG_3  = -DMANNING
# FLAG_4  = -DVESSEL
# FLAG_5  = -DMETEO
# FLAG_6  = -DWIND
# FLAG_7  = -DSEDIMENT
# FLAG_8  = -DCHECK_MASS_CONSERVATION
# FLAG_9  = -DTMP
# FLAG_10 = -DTRACKING
ifndef SEQU
  FLAG_PARALLEL = -DPARALLEL
endif
ifdef DEBUG
  FLAG_DEBUG = -DDEBUG
endif

EXEC = funwave

include Makefile-GNU.inc
#include Makefile-Intel.inc

CPPARGS = $(CPPFLAGS) $(FLAG_PRECISION) $(FLAG_CARTESIAN) $(FLAG_PARALLEL)\
 $(FLAG_1) $(FLAG_2) $(FLAG_3) $(FLAG_4) $(FLAG_5) $(FLAG_6) $(FLAG_7)\
 $(FLAG_8) $(FLAG_9) $(FLAG_10) $(FLAG_DEBUG)

#----------------FUNWAVE souce code-------------------
OBJS := main.o bc.o fluxes.o init.o io.o tridiagnal.o \
        breaker.o derivatives.o dispersion.o etauv_solver.o \
        sponge.o sources.o masks.o parallel.o statistics.o \
        wavemaker.o mixing.o nesting.o misc.o samples.o \
        mod_param.o mod_global.o mod_input.o mod_vessel.o mod_bathy_correction.o \
        mod_meteo.o mod_parallel_field_io.o mod_sediment.o mod_tide.o\
        mod_tracer.o

.PHONY: all

all: $(EXEC)

# link and generate the exec file
$(EXEC): $(OBJS)
	$(FC) $(DEBFLGS) $(INCS) -o $(EXEC) $(OBJS) $(LIBS)

# compile objective files
%.o : %.f90
	$(FC) $(FFLAGS) $(DEBFLGS) $(INCS) -o $@ -c $<

# compile objective to module has no action
%.mod : %.o

# pre-processors
%.f90 : $(SRC)/%.F
	$(CPP) $(CPPARGS) $< $@

mod_param.o: mod_param.f90
param.mod: mod_param.o
mod_global.o: mod_global.f90 param.mod
global.mod: mod_global.o
mod_input.o: mod_input.f90 param.mod global.mod
input_read.mod: mod_input.o
mod_vessel.o: mod_vessel.f90 param.mod global.mod input_read.mod
vessel_module.mod: mod_vessel.o
mod_bathy_correction.o: mod_bathy_correction.f90 param.mod global.mod input_read.mod
bathy_correction_module.mod: mod_bathy_correction.o
mod_meteo.o: mod_meteo.f90 param.mod global.mod input_read.mod
meteo_module.mod: mod_meteo.o
mod_parallel_field_io.o: mod_parallel_field_io.f90 global.mod
parallel_field_io.mod: mod_parallel_field_io.o
mod_sediment.o: mod_sediment.f90 param.mod global.mod input_read.mod
sediment_module.mod: mod_sediment.o
mod_tide.o: mod_tide.f90 param.mod global.mod input_read.mod
tide_module.mod: mod_tide.o
mod_tracer.o: mod_tracer.f90 param.mod global.mod input_read.mod
tracer.mod: mod_tracer.o

main.o: main.f90 global.mod vessel_module.mod sediment_module.mod meteo_module.mod tracer.mod tide_module.mod
bc.o: bc.f90 param.mod global.mod
fluxes.o: fluxes.f90 param.mod global.mod
init.o: init.f90 global.mod input_read.mod bathy_correction_module.mod
io.o: io.f90 param.mod global.mod input_read.mod vessel_module.mod sediment_module.mod meteo_module.mod tracer.mod parallel_field_io.mod
tridiagnal.o: tridiagnal.f90 param.mod global.mod
breaker.o: breaker.f90 param.mod global.mod
derivatives.o: derivatives.f90 param.mod
dispersion.o: dispersion.f90 global.mod
etauv_solver.o: etauv_solver.f90 global.mod vessel_module.mod
sponge.o: sponge.f90 param.mod global.mod
sources.o: sources.f90 global.mod vessel_module.mod meteo_module.mod
masks.o: masks.f90 global.mod vessel_module.mod
parallel.o: parallel.f90 param.mod global.mod
statistics.o: statistics.f90 global.mod
wavemaker.o: wavemaker.f90 param.mod global.mod
mixing.o: mixing.f90 param.mod global.mod
nesting.o: nesting.f90 global.mod input_read.mod
misc.o: misc.f90 param.mod global.mod meteo_module.mod
samples.o: samples.f90 param.mod global.mod

clean:
	rm -f *.o
	rm -f *.f90
	rm -f *.mod
	rm -f $(EXEC)


