# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../../..
SUBDIR=linsol/spfgmr/serial
C_SUBDIR=sunlinsol/spfgmr/serial

EXAMPLES = test_sunlinsol_spfgmr_serial.ml

include ../../linsol.mk

SPFGMR_SIZE ?= 100
SPFGMR_GSTYPE ?= 1
SPFGMR_MAXL ?= 500
SPFGMR_TOL ?= 0.0001
$(eval $(call EXECUTION_RULE,test_sunlinsol_spfgmr_serial,$$< $(SPFGMR_SIZE) $(SPFGMR_GSTYPE) $(SPFGMR_MAXL) $(SPFGMR_TOL) 0))

