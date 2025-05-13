# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../../..
SUBDIR=linsol/sptfqmr/parallel
C_SUBDIR=sunlinsol/sptfqmr/parallel

MPI_EXAMPLES = test_sunlinsol_sptfqmr_parallel.ml

include ../../linsol.mk

SPTFQMR_SIZE ?= 25
SPTFQMR_PRETYPE ?= 1
SPTFQMR_MAXL ?= 500
SPTFQMR_TOL ?= 0.0001
$(eval $(call EXECUTION_RULE,test_sunlinsol_sptfqmr_parallel,$(MPIRUN) -np 4 $$< $(SPTFQMR_SIZE) $(SPTFQMR_PRETYPE) $(SPTFQMR_MAXL) $(SPTFQMR_TOL) 0))

