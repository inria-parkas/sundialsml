# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../../..
SUBDIR=linsol/pcg/parallel
C_SUBDIR=sunlinsol/pcg/parallel

MPI_EXAMPLES = # $(if $(AT_LEAST_6_2_0),test_sunlinsol_pcg_parallel.ml)

include ../../linsol.mk

PCG_SIZE ?= 25
PCG_MAXL ?= 500
PCG_TOL ?= 0.0001
$(eval $(call EXECUTION_RULE,test_sunlinsol_pcg_parallel,$(MPIRUN) -np 4 $$< $(PCG_SIZE) $(PCG_MAXL) $(PCG_TOL) 0))

