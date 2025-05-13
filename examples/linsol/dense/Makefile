# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=linsol/dense
C_SUBDIR=sunlinsol/dense

EXAMPLES = test_sunlinsol_dense.ml

include ../linsol.mk

MATRIX_COLS ?= 500
$(eval $(call EXECUTION_RULE,test_sunlinsol_dense,$$< $(MATRIX_COLS) 0))

