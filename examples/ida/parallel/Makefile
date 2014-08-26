# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of IDAS.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless the
# latter are recompiled to link against IDAS.
USELIB ?= sundials_nosensi
SRCROOT=../../..
SUBDIR=ida/parallel
MPI_EXAMPLES = idaFoodWeb_kry_p.ml idaFoodWeb_kry_bbd_p.ml \
	       idaHeat2D_kry_p.ml idaHeat2D_kry_bbd_p.ml

include ../../examples.mk

MPIRUN = mpirun

# Tests with MPI.

idaFoodWeb_kry_p.byte.out: idaFoodWeb_kry_p.byte
	$(MPIRUN) -np 4 $< > $@
idaFoodWeb_kry_p.opt.out: idaFoodWeb_kry_p.opt
	$(MPIRUN) -np 4 $< > $@
idaFoodWeb_kry_p.sundials.out: $(EXAMPLESROOT)/ida/parallel/idaFoodWeb_kry_p
	$(MPIRUN) -np 4 $< > $@

idaFoodWeb_kry_bbd_p.byte.out: idaFoodWeb_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
idaFoodWeb_kry_bbd_p.opt.out: idaFoodWeb_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
idaFoodWeb_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/ida/parallel/idaFoodWeb_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

idaHeat2D_kry_bbd_p.byte.out: idaHeat2D_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
idaHeat2D_kry_bbd_p.opt.out: idaHeat2D_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
idaHeat2D_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/ida/parallel/idaHeat2D_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

idaHeat2D_kry_p.byte.out: idaHeat2D_kry_p.byte
	$(MPIRUN) -np 4 $< > $@
idaHeat2D_kry_p.opt.out: idaHeat2D_kry_p.opt
	$(MPIRUN) -np 4 $< > $@
idaHeat2D_kry_p.sundials.out: $(EXAMPLESROOT)/ida/parallel/idaHeat2D_kry_p
	$(MPIRUN) -np 4 $< > $@
