# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=linsol/lapackdense
C_SUBDIR=sunlinsol/lapackdense

LAPACK_EXAMPLES = test_sunlinsol_lapackdense.ml

include ../linsol.mk

MATRIX_COLS ?= 500
$(eval $(call EXECUTION_RULE,test_sunlinsol_lapackdense,$$< $(MATRIX_COLS) 0 0))

