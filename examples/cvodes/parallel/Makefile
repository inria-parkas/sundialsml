SRCROOT = ../../..
SUBDIR = cvodes/parallel
MPI_EXAMPLES = cvsAdvDiff_ASAp_non_p.ml		\
	       cvsAdvDiff_FSA_non_p.ml		\
	       cvsDiurnal_FSA_kry_p.ml		\
	       cvsAtmDisp_ASAi_kry_bbd_p.ml

include ../../examples.mk

# Tests with MPI and command-line arguments (make same choice as Sundials distrib)

$(eval $(call EXECUTION_RULE,cvsAdvDiff_ASAp_non_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,cvsAdvDiff_FSA_non_p,\
	$(MPIRUN) -np 4 $$< -sensi sim t))
$(eval $(call EXECUTION_RULE,cvsDiurnal_FSA_kry_p,\
	$(MPIRUN) -np 4 $$< -sensi sim t))
$(eval $(call EXECUTION_RULE,cvsAtmDisp_ASAi_kry_bbd_p,\
	$(MPIRUN) -np 4 $$< -sensi sim t))
