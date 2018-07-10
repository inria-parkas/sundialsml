# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=matrix/arraydense
C_SUBDIR=sunmatrix/dense

# A bit hackish: arraydense should be compared to dense, but the
# latter lives in a separate directory, making it hard to accurately
# track the dependency.  We reproduce all the rules here for building
# dense, but disable dense so that only the parts needed for building
# arraydense are live.
ENABLED_EXAMPLES = test_sunmatrix_arraydense.ml
EXAMPLES         = test_sunmatrix_arraydense.ml test_sunmatrix_dense.ml

include ../matrix.mk

MATRIX_ROWS ?= 1000
MATRIX_COLS ?= 500
$(eval $(call EXECUTION_RULE,test_sunmatrix_arraydense,$$< $(MATRIX_ROWS) $(MATRIX_COLS) 0,custom-sundials-rule))
$(eval $(call EXECUTION_RULE,test_sunmatrix_dense,$$< $(MATRIX_ROWS) $(MATRIX_COLS) 0))

test_sunmatrix_arraydense.sundials.out: test_sunmatrix_dense.sundials.out
	sed 's/^Dense matrix test:/ArrayDense matrix test:/' $< > $@
