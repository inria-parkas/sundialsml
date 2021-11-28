# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=nvector/mpiplusx

MPI_EXAMPLES = \
	$(if $(AT_LEAST_5_0_0), test_nvector_mpiplusx.ml)

include ../nvector.mk

NVECTOR_SIZE ?= 5000
$(eval $(call EXECUTION_RULE,test_nvector_mpiplusx, \
	 $(MPIRUN) -np 4 $(UTILS)/unmute.sh 0 $$< $(NVECTOR_SIZE) 0))

