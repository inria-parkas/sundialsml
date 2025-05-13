# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=linsol/superlumt
C_SUBDIR=sunlinsol/superlumt

EXAMPLES = $(if $(SUPERLUMT_ENABLED),test_sunlinsol_superlumt.ml)

include ../linsol.mk

MATRIX_SIZE ?= 500
MATRIX_TYPE ?= 1
NUMTHREADS ?= 1
$(eval $(call EXECUTION_RULE,test_sunlinsol_superlumt,$$< $(MATRIX_SIZE) $(MATRIX_TYPE) $(NUMTHREADS) 0))

