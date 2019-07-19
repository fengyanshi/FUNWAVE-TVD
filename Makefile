#-----------BEGIN MAKEFILE---------------------------
FUNWAVE_DIR = .
WORK_DIR    = funwave-work
COMPILER    = gnu
PARALLEL    = true
EXEC        = funwave
PRECISION   = single

#-----------DEFINE FLAGS-----------------------------
#         uncomment to choose the model
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

#----------------uncommon options---------------------
DEF_FC      = 
DEF_FC_FLAG = 
SPHERICAL   = false
MPI         = openmpi
DEBUG       = false
INCS        = $(IOINCS) $(GOTMINCS)
LIBS        = $(PV3LIB) $(CLIB)  $(PARLIB) $(IOLIBS) $(MPILIB) $(GOTMLIB)
CLIB        =
MDEPFLAGS   = --cpp --fext=f90 --file=-
RANLIB      = ranlib

#----------include the essential makefiles------------
include $(FUNWAVE_DIR)/GNUMake/Essential/Make_Essential

##-----------------------------------------------------
##      Instructions for Makefile
##-----------------------------------------------------
 
##--------Make options (provided by Make_Essential)----
# make:  
#      Create the $(WORK_DIR) directory and the exectutable
# make clean:
#      Clean the "build" directory in "$(WORK_DIR)" directory
# make clean-exe:
#      Clean the "build" directory in "$(WORK_DIR)" directory and the executable
# make clobber:
#      Clean up the whole $(WORK_DIR) directory
# make chec-env:
#      Print the compier version and mpi version
# make print-foo
#      Check the value of "foo" in Makefile (Or Makefile_Essential)
#      For example, "make print $(EXEC)" and you will see the final $(EXEC) name

##---------Notes to the Makefils variables------------
# FUNWAVE_DIR: 
#      The path to FUNWAVE directory, can be either absolute path or reference path. 
#      To set to current directory, set      FUNWAVE_DIR = .
# WORK_DIR:
#      The path to the work directory. A new directory will be created named $(WORK_DIR)
#      To set to current directory, set      FUNWAVE_DIR = .
# COMPILER:
#      Support list: "gnu", "intel"
# PARALLEL:
#      eithrt "true" or "false"
# EXEC: 
#      The name for exectutable. By default, the name of the executable will be "funwave+$(SUFFIX)".
#      SUFFIX is a self-explainary string depending on the options above
# PRESICION:
#      either "single" or "double"
# 
# DEF_FC:
#      Use Defined compiler name for specific machine.
#      Left for empty for dafault.
# DEF_FC_FLAG:
#      User Defined comiler flags for $(DEF_FFC)
#      Left for empty for dafault.
# SPHIRICAL:
#      "true" or "false". set to "true" to enable SPHERICAL coordinate.
# MPI:
#      Depend on the mpi version. In the settings of Make_Essential,
#      if PARALLEL=true, MPI=openmpi or MPImpich, FC = mpif90 (For any COMPILER)
#      if PARALLEL=true, COMPILER=intel, MPI=intelmpi, FC = mpiifort
# DEBUG:
#      either "true" or "false"
#      If "true" uses some debug flags, for example -Wall for gnu compiler
#      If "false", no debug flags used. Use optimization flag (-O2 for intel, -O3 for gnu)

## To use like the 3.3 or earlier version of FUNWAVE-TVD
# Go to "src" directory and set and set both "FUNWAVE_DIR" and "WORK_DIR" as ".".
# Here is an example for vessel case (uncomment all below). and you will see
#-----------BEGIN MAKEFILE---------------------------
# FUNWAVE_DIR = .
# WORK_DIR    = .
# COMPILER    = gnu
# PARALLEL    = true
# EXEC        = funwave
# PRECISION   = double

# #-----------DEFINE FLAGS-----------------------------
# #         uncomment to choose the model
# # FLAG_1  = -DCOUPLING
# # FLAG_2  = -DZALPHA
# # FLAG_3  = -DMANNING
# # FLAG_4  = -DVESSEL
# # FLAG_5  = -DMETEO
# # FLAG_6  = -DWIND
# # FLAG_7  = -DSEDIMENT
# FLAG_8  = -DCHECK_MASS_CONSERVATION
# # FLAG_9  = -DTMP
# # FLAG_10 = -DTRACKING

# #----------------uncommon options---------------------
# SPHERICAL = false
# MPI       = openmpi
# DEBUG     = false
# MDEPFLAGS = --cpp --fext=f90 --file=-
# RANLIB = ranlib
