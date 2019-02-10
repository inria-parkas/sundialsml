# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=matrix/sparse
C_SUBDIR=sunmatrix/sparse

EXAMPLES = test_sunmatrix_sparse.ml

include ../matrix.mk

MATRIX_ROWS ?= 500
MATRIX_COLS ?= 500
MATRIX_TYPE ?= 0
$(eval $(call EXECUTION_RULE,test_sunmatrix_sparse,$$< $(MATRIX_ROWS) $(MATRIX_COLS) $(MATRIX_TYPE) 0))

