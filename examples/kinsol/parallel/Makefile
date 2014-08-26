# Invoke `make' with USELIB=sundials to run the tests with the
# sundials.cma rather than sundials_nosensi.cma.  It should make no
# difference to the test results.
USELIB ?= sundials_nosensi
SRCROOT = ../../..
SUBDIR = kinsol/parallel

MPI_EXAMPLES = kinFoodWeb_kry_p.ml		\
	       kinFoodWeb_kry_bbd_p.ml

include ../../examples.mk

MPIRUN = mpirun

# Tests with MPI.

kinFoodWeb_kry_p.byte.out: kinFoodWeb_kry_p.byte
	$(MPIRUN) -np 4 $< > $@
kinFoodWeb_kry_p.opt.out: kinFoodWeb_kry_p.opt
	$(MPIRUN) -np 4 $< > $@
kinFoodWeb_kry_p.sundials.out: $(EXAMPLESROOT)/kinsol/parallel/kinFoodWeb_kry_p
	$(MPIRUN) -np 4 $< > $@

kinFoodWeb_kry_bbd_p.byte.out: kinFoodWeb_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
kinFoodWeb_kry_bbd_p.opt.out: kinFoodWeb_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
kinFoodWeb_kry_bbd_p.sundials.out: \
    $(EXAMPLESROOT)/kinsol/parallel/kinFoodWeb_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@
