# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=arkode/C_manyvector

EXAMPLES_SINCE_5_0_0 = ark_brusselator1D_manyvec.ml

EXAMPLES = $(if $(AT_LEAST_5_0_0),$(EXAMPLES_SINCE_5_0_0))

FILES_TO_CLEAN = bruss_mesh.txt bruss_u.txt bruss_v.txt bruss_w.txt

include ../../examples.mk
