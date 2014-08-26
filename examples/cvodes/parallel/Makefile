SRCROOT = ../../..
SUBDIR = cvodes/parallel
MPI_EXAMPLES = cvsAdvDiff_ASAp_non_p.ml		\
	       cvsAdvDiff_FSA_non_p.ml		\
	       cvsDiurnal_FSA_kry_p.ml		\
	       cvsAtmDisp_ASAi_kry_bbd_p.ml

include ../../examples.mk

MPIRUN = mpirun

# Tests with command-line arguments (make same choice as Sundials distrib)

cvsAdvDiff_ASAp_non_p.byte.out: cvsAdvDiff_ASAp_non_p.byte
	$(MPIRUN) -np 4 $< > $@
cvsAdvDiff_ASAp_non_p.opt.out: cvsAdvDiff_ASAp_non_p.opt
	$(MPIRUN) -np 4 $< > $@
cvsAdvDiff_ASAp_non_p.sundials.out: $(EXAMPLESROOT)/cvodes/parallel/cvsAdvDiff_ASAp_non_p
	$(MPIRUN) -np 4 $< > $@

cvsAdvDiff_FSA_non_p.byte.out: cvsAdvDiff_FSA_non_p.byte
	$(MPIRUN) -np 4 $< -sensi sim t > $@
cvsAdvDiff_FSA_non_p.opt.out: cvsAdvDiff_FSA_non_p.opt
	$(MPIRUN) -np 4 $< -sensi sim t > $@
cvsAdvDiff_FSA_non_p.sundials.out: $(EXAMPLESROOT)/cvodes/parallel/cvsAdvDiff_FSA_non_p
	$(MPIRUN) -np 4 $< -sensi sim t > $@

cvsDiurnal_FSA_kry_p.byte.out: cvsDiurnal_FSA_kry_p.byte
	$(MPIRUN) -np 4 $< -sensi sim t > $@
cvsDiurnal_FSA_kry_p.opt.out: cvsDiurnal_FSA_kry_p.opt
	$(MPIRUN) -np 4 $< -sensi sim t > $@
cvsDiurnal_FSA_kry_p.sundials.out: $(EXAMPLESROOT)/cvodes/parallel/cvsDiurnal_FSA_kry_p
	$(MPIRUN) -np 4 $< -sensi sim t > $@

cvsAtmDisp_ASAi_kry_bbd_p.byte.out: cvsAtmDisp_ASAi_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< -sensi sim t > $@
cvsAtmDisp_ASAi_kry_bbd_p.opt.out: cvsAtmDisp_ASAi_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< -sensi sim t > $@
cvsAtmDisp_ASAi_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/cvodes/parallel/cvsAtmDisp_ASAi_kry_bbd_p
	$(MPIRUN) -np 4 $< -sensi sim t > $@
