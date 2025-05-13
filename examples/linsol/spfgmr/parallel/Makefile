# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../../..
SUBDIR=linsol/spfgmr/parallel
C_SUBDIR=sunlinsol/spfgmr/parallel

MPI_EXAMPLES = test_sunlinsol_spfgmr_parallel.ml

include ../../linsol.mk

SPFGMR_SIZE ?= 25
SPFGMR_GSTYPE ?= 1
SPFGMR_MAXL ?= 500
SPFGMR_TOL ?= 0.0001
$(eval $(call EXECUTION_RULE,test_sunlinsol_spfgmr_parallel,$(MPIRUN) -np 4 $$< $(SPFGMR_SIZE) $(SPFGMR_GSTYPE) $(SPFGMR_MAXL) $(SPFGMR_TOL) 0))

