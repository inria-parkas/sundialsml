# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=linsol/klu
C_SUBDIR=sunlinsol/klu

EXAMPLES = $(if $(KLU_ENABLED),test_sunlinsol_klu.ml)

include ../linsol.mk

MATRIX_COLS ?= 500
MATRIX_TYPE ?= 0
$(eval $(call EXECUTION_RULE,test_sunlinsol_klu,$$< $(MATRIX_COLS) $(MATRIX_TYPE) 0))

