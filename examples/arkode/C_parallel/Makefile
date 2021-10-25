# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=arkode/C_parallel

EXAMPLES_SINCE_5_4_0 = ark_brusselator1D_task_local_nls.ml

MPI_EXAMPLES = ark_diurnal_kry_p.ml \
	       ark_diurnal_kry_bbd_p.ml \
	   	$(if $(AT_LEAST_5_4_0),$(EXAMPLES_SINCE_5_4_0))

include ../../examples.mk

$(eval $(call EXECUTION_RULE,ark_diurnal_kry_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,ark_diurnal_kry_bbd_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,ark_brusselator1D_task_local_nls,$(MPIRUN) -np 4 $$<))

