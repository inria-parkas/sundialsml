# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=nvector/serial

EXAMPLES = test_nvector_serial.ml

include ../nvector.mk

NVECTOR_SIZE ?= 100000
$(eval $(call EXECUTION_RULE,test_nvector_serial,$$< $(NVECTOR_SIZE) 0))
