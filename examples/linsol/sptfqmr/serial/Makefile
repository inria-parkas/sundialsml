# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../../..
SUBDIR=linsol/sptfqmr/serial
C_SUBDIR=sunlinsol/sptfqmr/serial

EXAMPLES = test_sunlinsol_sptfqmr_serial.ml

include ../../linsol.mk

SPTFQMR_SIZE ?= 100
SPTFQMR_PRETYPE ?= 1
SPTFQMR_MAXL ?= 500
SPTFQMR_TOL ?= 0.0001
$(eval $(call EXECUTION_RULE,test_sunlinsol_sptfqmr_serial,$$< $(SPTFQMR_SIZE) $(SPTFQMR_PRETYPE) $(SPTFQMR_MAXL) $(SPTFQMR_TOL) 0))

