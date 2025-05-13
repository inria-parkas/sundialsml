# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../../..
SUBDIR=linsol/spgmr/parallel
C_SUBDIR=sunlinsol/spgmr/parallel

MPI_EXAMPLES = test_sunlinsol_spgmr_parallel.ml

include ../../linsol.mk

SPGMR_SIZE ?= 25
SPGMR_GSTYPE ?= 1
SPGMR_PRETYPE ?= 1
SPGMR_MAXL ?= 500
SPGMR_TOL ?= 0.0001
$(eval $(call EXECUTION_RULE,test_sunlinsol_spgmr_parallel,$(MPIRUN) -np 4 $$< $(SPGMR_SIZE) $(SPGMR_GSTYPE) $(SPGMR_PRETYPE) $(SPGMR_MAXL) $(SPGMR_TOL) 0))

