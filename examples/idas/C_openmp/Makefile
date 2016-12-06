# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of IDAS.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless the
# latter are recompiled to link against IDAS.
USELIB ?= sundials
SRCROOT=../../..
SUBDIR=idas/C_openmp

OPENMP_EXAMPLES = $(if $(AT_LEAST_2_7),idasFoodWeb_bnd_omp.ml) \
		  $(if $(AT_LEAST_2_7),idasFoodWeb_kry_omp.ml)

include ../../examples.mk

# Run with 8 threads.
$(eval $(call EXECUTION_RULE,idasFoodWeb_bnd_omp,$$< 8))
# Run with 1 thread since the Sundials version is not deterministic otherwise.
$(eval $(call EXECUTION_RULE,idasFoodWeb_kry_omp,$$< 1))
