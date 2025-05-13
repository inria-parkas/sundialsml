# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../../..
SUBDIR=linsol/spbcgs/serial
C_SUBDIR=sunlinsol/spbcgs/serial

EXAMPLES = test_sunlinsol_spbcgs_serial.ml

include ../../linsol.mk

SPBCGS_SIZE ?= 100
SPBCGS_PRETYPE ?= 1
SPBCGS_MAXL ?= 500
SPBCGS_TOL ?= 0.0001
$(eval $(call EXECUTION_RULE,test_sunlinsol_spbcgs_serial,$$< $(SPBCGS_SIZE) $(SPBCGS_PRETYPE) $(SPBCGS_MAXL) $(SPBCGS_TOL) 0))

