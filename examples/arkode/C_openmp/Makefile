# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=arkode/C_openmp

OPENMP_EXAMPLES = ark_brusselator1D_omp.ml

FILES_TO_CLEAN = bruss_mesh.txt bruss_u.txt bruss_v.txt bruss_w.txt \
		 heat1D.txt heat_mesh.txt solution.txt

include ../../examples.mk

# Run with 4 threads.
$(eval $(call EXECUTION_RULE,ark_brusselator1D_omp,$$< 4))

