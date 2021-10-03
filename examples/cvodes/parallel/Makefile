SRCROOT = ../../..
SUBDIR = cvodes/parallel
MPI_EXAMPLES = cvsAdvDiff_ASAp_non_p.ml		\
	       cvsAdvDiff_FSA_non_p.ml		\
	       cvsDiurnal_FSA_kry_p.ml		\
	       cvsAtmDisp_ASAi_kry_bbd_p.ml

FILES_TO_CLEAN = grad.m grad000.m grad001.m grad002.m grad003.m

include ../../examples.mk

# Tests with MPI and command-line arguments (make same choice as Sundials distrib)

$(eval $(call EXECUTION_RULE,cvsAdvDiff_ASAp_non_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,cvsAdvDiff_FSA_non_p,\
	$(MPIRUN) -np 4 $$< -sensi sim t))
$(eval $(call EXECUTION_RULE,cvsDiurnal_FSA_kry_p,\
	$(MPIRUN) -np 4 $$< -sensi sim t))
ifeq ($(AT_LEAST_2_7),)
    $(eval $(call EXECUTION_RULE,cvsAtmDisp_ASAi_kry_bbd_p,\
	    $(MPIRUN) -np 4 $$< -sensi sim t))
else
ifeq ($(AT_LEAST_5_0_0),)
    $(eval $(call EXECUTION_RULE,cvsAtmDisp_ASAi_kry_bbd_p,\
	    $(MPIRUN) -np 8 $$< -sensi sim t))
else
    $(eval $(call EXECUTION_RULE,cvsAtmDisp_ASAi_kry_bbd_p,\
	    $(MPIRUN) -np 4 $$< -sensi sim t))
endif
endif
