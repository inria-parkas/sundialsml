# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=arkode/C_openmp

OPENMP_EXAMPLES_SINCE_2_7_0 = ark_brusselator1D_omp.ml

OPENMP_EXAMPLES_SINCE_4_0_0 = ark_heat1D_omp.ml

OPENMP_EXAMPLES = $(OPENMP_EXAMPLES_SINCE_2_7_0) \
		  $(if $(AT_LEAST_4_0_0),$(OPENMP_EXAMPLES_SINCE_4_0_0))

FILES_TO_CLEAN = bruss_mesh.txt bruss_u.txt bruss_v.txt bruss_w.txt \
		 heat1D.txt heat_mesh.txt solution.txt

include ../../examples.mk

# Run with 4 threads.
$(eval $(call EXECUTION_RULE,ark_brusselator1D_omp,$$< 4))

