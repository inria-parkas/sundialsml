# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=arkode/C_parallel

EXAMPLES = ark_diurnal_kry_p.ml

include ../../examples.mk

$(eval $(call EXECUTION_RULE,ark_diurnal_kry_p,$(MPIRUN) -np 4 $$<))
