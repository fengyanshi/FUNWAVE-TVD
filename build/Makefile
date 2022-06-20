SRC=../src
DST=../build
# FLAG_PRECISION = -DDOUBLE_PRECISION 
FLAG_CARTESIAN = -DCARTESIAN
FLAG_1  = -DCOUPLING
FLAG_2  = -DZALPHA
FLAG_3  = -DMANNING
FLAG_4  = -DVESSEL
FLAG_5  = -DMETEO
FLAG_6  = -DWIND
FLAG_7  = -DSEDIMENT
FLAG_8  = -DCHECK_MASS_CONSERVATION
FLAG_9  = -DTMP
FLAG_10 = -DTRACKING
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

# pre-processors
%.f90 : $(SRC)/%.F
	$(CPP) $(CPPARGS) $< $@

#mod_param.f90: $(SRC)/mod_param.F
#mod_global.f90: $(SRC)/mod_global.F
#mod_input.f90: $(SRC)/mod_input.F
#mod_vessel.f90: $(SRC)/mod_vessel.F
#mod_bathy_correction.f90: $(SRC)/mod_bathy_correction.F
#mod_meteo.f90: $(SRC)/mod_meteo.F
#mod_parallel_field_io.f90: $(SRC)/mod_parallel_field_io.F
#mod_sediment.f90: $(SRC)/mod_sediment.F
#mod_tide.f90: $(SRC)/mod_tide.F
#mod_tracer.f90: $(SRC)/mod_tracer.F

#main.f90: $(SRC)/main.F
#bc.f90: $(SRC)/bc.F
#fluxes.f90: $(SRC)/fluxes.F
#init.f90: $(SRC)/init.F
#io.f90: $(SRC)/io.F
#tridiagnal.f90: $(SRC)/tridiagnal.F
#breaker.f90: $(SRC)/breaker.F
#derivatives.f90: $(SRC)/derivatives.F
#dispersion.f90: $(SRC)/dispersion.F
#etauv_solver.f90: $(SRC)/etauv_solver.F
#sponge.f90: $(SRC)/sponge.F
#sources.f90: $(SRC)/sources.F
#masks.f90: $(SRC)/masks.F
#parallel.f90: $(SRC)/parallel.F
#statistics.f90: $(SRC)/statistics.F
#wavemaker.f90: $(SRC)/wavemaker.F
#mixing.f90: $(SRC)/mixing.F
#nesting.f90: $(SRC)/nesting.F
#misc.f90: $(SRC)/misc.F
#samples.f90: $(SRC)/samples.F

mod_param.o: mod_param.f90
mod_global.o: mod_global.f90 mod_param.o
mod_input.o: mod_input.f90 mod_param.o mod_global.o
mod_vessel.o: mod_vessel.f90 mod_param.o mod_global.o mod_input.o
mod_bathy_correction.o: mod_bathy_correction.f90 mod_param.o mod_global.o mod_input.o
mod_meteo.o: mod_meteo.f90 mod_param.o mod_global.o mod_input.o
mod_parallel_field_io.o: mod_parallel_field_io.f90 mod_global.o
mod_sediment.o: mod_sediment.f90 mod_param.o mod_global.o mod_input.o
mod_tide.o: mod_tide.f90 mod_param.o mod_global.o mod_input.o
mod_tracer.o: mod_tracer.f90 mod_param.o mod_global.o mod_input.o

main.o: main.f90 mod_global.o mod_vessel.o mod_sediment.o mod_meteo.o mod_tracer.o mod_tide.o
bc.o: bc.f90 mod_param.o mod_global.o
fluxes.o: fluxes.f90 mod_param.o mod_global.o
init.o: init.f90 mod_global.o mod_input.o mod_bathy_correction.o
io.o: io.f90 mod_param.o mod_global.o mod_input.o mod_vessel.o mod_sediment.o mod_meteo.o mod_tracer.o mod_parallel_field_io.o
tridiagnal.o: tridiagnal.f90 mod_param.o mod_global.o
breaker.o: breaker.f90 mod_param.o mod_global.o
derivatives.o: derivatives.f90 mod_param.o
dispersion.o: dispersion.f90 mod_global.o
etauv_solver.o: etauv_solver.f90 mod_global.o mod_vessel.o
sponge.o: sponge.f90 mod_param.o mod_global.o
sources.o: sources.f90 mod_global.o mod_vessel.o mod_meteo.o
masks.o: masks.f90 mod_global.o mod_vessel.o
parallel.o: parallel.f90 mod_param.o mod_global.o
statistics.o: statistics.f90 mod_global.o
wavemaker.o: wavemaker.f90 mod_param.o mod_global.o
mixing.o: mixing.f90 mod_param.o mod_global.o
nesting.o: nesting.f90 mod_global.o mod_input.o
misc.o: misc.f90 mod_param.o mod_global.o mod_meteo.o
samples.o: samples.f90 mod_param.o mod_global.o

clean:
	rm -f *.o
	rm -f *.f90
	rm -f *.mod
	rm -f $(EXEC)


