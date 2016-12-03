# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=cvode/C_openmp

OPENMP_EXAMPLES = $(if $(AT_LEAST_2_7),cvAdvDiff_bnd_omp.ml)

include ../../examples.mk

# Run with 4 threads.
$(eval $(call EXECUTION_RULE,cvAdvDiff_bnd_omp,$$< 4))
