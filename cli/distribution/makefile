#==================================================
# Makefile to build the 'distribution' CLI utility
#==================================================
#

# Define the Fortran Compiler and its options
FC  = gfortran
FLAGS = -fcheck=all -Wall -O3 # DEBUG: -fcheck=all -Wall -g -O0

# Define directories
DMSL_DIR = ../../../miniDMSL/src/
BMSL_DIR = ../../src/

# Define DMSL files
DMSL = \
	$(DMSL_DIR)kinds_dmsl_kit.f90              \
	$(DMSL_DIR)types_dmsl_kit.f90              \
	$(DMSL_DIR)utilities_dmsl_kit.f90          \
	$(DMSL_DIR)mt19937-64.f95                  \
	$(DMSL_DIR)uniran1_minidmsl_mod.f90        \
	$(DMSL_DIR)linalg_dmsl_kit.f90             \
	$(DMSL_DIR)numerix_dmsl_kit.f90

# Define BMSL files
BMSL = \
	$(BMSL_DIR)EmpiricalStats_tools.f90          \
	$(BMSL_DIR)Distribution_tools.f90 

# Define distribution files
dist = main_distribution.f90 
	
all:build_exe

build_exe: compile_all
	$(FC)  $(FLAGS) $(DMSL) $(BMSL) $(dist) -o distribution

compile_all: compile_bmsl
	$(FC) $(FLAGS) -c $(dist)

compile_bmsl: compile_dmsl
	$(FC) $(FLAGS) -c $(BMSL)

compile_dmsl:
	$(FC) $(FLAGS) -c $(DMSL)





