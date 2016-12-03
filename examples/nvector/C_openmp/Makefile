# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..

SUBDIR=nvector/C_openmp
ifeq ($(AT_LEAST_2_7),)
    C_SUBDIR=nvector/C_openmp
else
    C_SUBDIR=nvector/openmp
endif

OPENMP_EXAMPLES = test_nvector_openmp.ml

include ../nvector.mk

NVECTOR_SIZE ?= 100000
NUM_THREADS ?= 4
$(eval $(call EXECUTION_RULE,test_nvector_openmp,\
	$$< $(NVECTOR_SIZE) $(NUM_THREADS) 0))
