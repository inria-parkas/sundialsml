# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=cvode/C_mpimanyvector

EXAMPLES_SINCE_5_0_0 = cvDiurnal_kry_mpimanyvec.ml

MPI_EXAMPLES = $(if $(AT_LEAST_5_0_0),$(EXAMPLES_SINCE_5_0_0))

include ../../examples.mk

$(eval $(call EXECUTION_RULE,cvDiurnal_kry_mpimanyvec,$(MPIRUN) -np 4 $$<))
