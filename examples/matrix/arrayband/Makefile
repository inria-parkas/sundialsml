# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=matrix/arrayband
C_SUBDIR=sunmatrix/band

# A bit hackish: arrayband should be compared to band, but the latter
# lives in a separate directory, making it hard to accurately track
# the dependency.  We reproduce all the rules here for building band,
# but disable band so that only the parts needed for building
# arrayband are live.
ENABLED_EXAMPLES = test_sunmatrix_arrayband.ml
EXAMPLES         = test_sunmatrix_arrayband.ml test_sunmatrix_band.ml

include ../matrix.mk

MATRIX_COLS ?= 5000
MATRIX_UBAND ?= 100
MATRIX_LBAND ?= 100
$(eval $(call EXECUTION_RULE,test_sunmatrix_arrayband,$$< $(MATRIX_COLS) $(MATRIX_UBAND) $(MATRIX_LBAND) 0,custom-sundials-rule))
$(eval $(call EXECUTION_RULE,test_sunmatrix_band,$$< $(MATRIX_COLS) $(MATRIX_UBAND) $(MATRIX_LBAND) 0))

test_sunmatrix_arrayband.sundials.out: test_sunmatrix_band.sundials.out
	sed 's/^Band matrix test:/ArrayBand matrix test:/' $< > $@
