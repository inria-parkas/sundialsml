# Invoke `make' with USELIB=sundials to run the tests with the
# sundials.cma rather than sundials_nosensi.cma.  It should make no
# difference to the test results.
USELIB ?= sundials_nosensi
SRCROOT = ../../..
SUBDIR = kinsol/parallel

MPI_EXAMPLES = kinFoodWeb_kry_p.ml		\
	       kinFoodWeb_kry_bbd_p.ml

include ../../examples.mk

# Tests with MPI.

$(eval $(call EXECUTION_RULE,kinFoodWeb_kry_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,kinFoodWeb_kry_bbd_p,$(MPIRUN) -np 4 $$<))
