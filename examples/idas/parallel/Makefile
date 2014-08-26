SRCROOT = ../../..
SUBDIR = idas/parallel

MPI_EXAMPLES = idasHeat2D_FSA_kry_bbd_p.ml	\
	       idasBruss_kry_bbd_p.ml		\
	       idasBruss_FSA_kry_bbd_p.ml	\
	       idasBruss_ASAp_kry_bbd_p.ml

include ../../examples.mk

# Tests with MPI and command-line arguments (make same choice as Sundials distrib)

MPIRUN = mpirun

# This test cannot be run with -sensi sim.
idasHeat2D_FSA_kry_bbd_p.byte.out: idasHeat2D_FSA_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< -sensi stg t > $@
idasHeat2D_FSA_kry_bbd_p.opt.out: idasHeat2D_FSA_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< -sensi stg t > $@
idasHeat2D_FSA_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/idas/parallel/idasHeat2D_FSA_kry_bbd_p
	$(MPIRUN) -np 4 $< -sensi stg t > $@

idasBruss_kry_bbd_p.byte.out: idasBruss_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
idasBruss_kry_bbd_p.opt.out: idasBruss_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
idasBruss_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/idas/parallel/idasBruss_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

idasBruss_FSA_kry_bbd_p.byte.out: idasBruss_FSA_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
idasBruss_FSA_kry_bbd_p.opt.out: idasBruss_FSA_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
idasBruss_FSA_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/idas/parallel/idasBruss_FSA_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

idasBruss_ASAp_kry_bbd_p.byte.out: idasBruss_ASAp_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
idasBruss_ASAp_kry_bbd_p.opt.out: idasBruss_ASAp_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
idasBruss_ASAp_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/idas/parallel/idasBruss_ASAp_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@
