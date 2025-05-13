# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=linsol/band
C_SUBDIR=sunlinsol/band

EXAMPLES = test_sunlinsol_band.ml

include ../linsol.mk

MATRIX_COLS ?= 5000
MATRIX_UBAND ?= 100
MATRIX_LBAND ?= 100
$(eval $(call EXECUTION_RULE,test_sunlinsol_band,$$< $(MATRIX_COLS) $(MATRIX_UBAND) $(MATRIX_LBAND) 0))

