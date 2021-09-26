# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=nvector/mpimanyvector

MPI_EXAMPLES = \
	$(if $(AT_LEAST_5_0_0), test_nvector_mpimanyvector_serial.ml) \
	$(if $(AT_LEAST_5_0_0), test_nvector_mpimanyvector_parallel1.ml) \
	$(if $(AT_LEAST_5_0_0), test_nvector_mpimanyvector_parallel2.ml)

include ../nvector.mk

NVECTOR_SIZE1 ?= 50000
NVECTOR_SIZE2 ?= 20000
$(eval $(call EXECUTION_RULE,test_nvector_mpimanyvector_serial, \
	 $$< $(NVECTOR_SIZE1) $(NVECTOR_SIZE2) 0))
$(eval $(call EXECUTION_RULE,test_nvector_mpimanyvector_parallel1, \
	 $(MPIRUN) -np 4 $$< $(NVECTOR_SIZE1) $(NVECTOR_SIZE2) 0))
$(eval $(call EXECUTION_RULE,test_nvector_mpimanyvector_parallel2, \
	 $(MPIRUN) -np 4 $$< $(NVECTOR_SIZE1) $(NVECTOR_SIZE2) 0))
