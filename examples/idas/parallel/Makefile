SRCROOT = ../../..
SUBDIR = idas/parallel

MPI_EXAMPLES = idasHeat2D_FSA_kry_bbd_p.ml	\
	       idasBruss_kry_bbd_p.ml		\
	       idasBruss_FSA_kry_bbd_p.ml	\
	       idasBruss_ASAp_kry_bbd_p.ml

include ../../examples.mk

# Tests with MPI and command-line arguments (make same choice as Sundials distrib)

# This test cannot be run with -sensi sim.
$(eval $(call EXECUTION_RULE,idasHeat2D_FSA_kry_bbd_p,\
	$(MPIRUN) -np 4 $$< -sensi stg t))
$(eval $(call EXECUTION_RULE,idasBruss_kry_bbd_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,idasBruss_FSA_kry_bbd_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,idasBruss_ASAp_kry_bbd_p,$(MPIRUN) -np 4 $$<))
