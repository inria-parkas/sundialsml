# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_nosensi
SRCROOT=../../..
SUBDIR=cvode/parallel

MPI_EXAMPLES = cvDiurnal_kry_p.ml \
	       cvDiurnal_kry_bbd_p.ml \
	       cvAdvDiff_non_p.ml

include ../../examples.mk

MPIRUN = mpirun

cvDiurnal_kry_p.byte.out: cvDiurnal_kry_p.byte
	$(MPIRUN) -np 4 $< > $@
cvDiurnal_kry_p.opt.out: cvDiurnal_kry_p.opt
	$(MPIRUN) -np 4 $< > $@
cvDiurnal_kry_p.sundials.out: $(EXAMPLESROOT)/cvode/parallel/cvDiurnal_kry_p
	$(MPIRUN) -np 4 $< > $@

cvDiurnal_kry_bbd_p.byte.out: cvDiurnal_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
cvDiurnal_kry_bbd_p.opt.out: cvDiurnal_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
cvDiurnal_kry_bbd_p.sundials.out: \
			$(EXAMPLESROOT)/cvode/parallel/cvDiurnal_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

cvAdvDiff_non_p.byte.out: cvAdvDiff_non_p.byte
	$(MPIRUN) -np 6 $< > $@
cvAdvDiff_non_p.opt.out: cvAdvDiff_non_p.opt
	$(MPIRUN) -np 6 $< > $@
cvAdvDiff_non_p.sundials.out: \
			$(EXAMPLESROOT)/cvode/parallel/cvAdvDiff_non_p
	$(MPIRUN) -np 6 $< > $@
