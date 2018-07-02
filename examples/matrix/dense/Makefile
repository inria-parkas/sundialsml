# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=matrix/dense
C_SUBDIR=sunmatrix/dense

EXAMPLES = test_sunmatrix_dense.ml

include ../matrix.mk

MATRIX_ROWS ?= 1000
MATRIX_COLS ?= 500
$(eval $(call EXECUTION_RULE,test_sunmatrix_dense,$$< $(MATRIX_ROWS) $(MATRIX_COLS) 0))

