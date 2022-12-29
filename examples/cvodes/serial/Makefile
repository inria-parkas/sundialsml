SRCROOT=../../..
SUBDIR=cvodes/serial

EXAMPLES = cvsAdvDiff_ASAi_bnd.ml \
	   cvsAdvDiff_FSA_non.ml  \
	   cvsDiurnal_FSA_kry.ml  \
	   cvsFoodWeb_ASAi_kry.ml \
	   cvsFoodWeb_ASAp_kry.ml \
	   cvsHessian_ASA_FSA.ml  \
	   cvsRoberts_ASAi_dns.ml \
	   $(if $(OK_CVodeGetSensNumNonlinSolvIters),cvsRoberts_FSA_dns.ml) \
	   $(if $(and $(KLU_ENABLED),$(AT_LEAST_2_7)),cvsRoberts_ASAi_klu.ml) \
	   $(if $(and $(KLU_ENABLED),\
	   	      $(AT_LEAST_2_7),\
		      $(OK_CVodeGetSensNumNonlinSolvIters)),cvsRoberts_FSA_klu.ml) \
	   $(if $(and $(SUPERLUMT_ENABLED),$(AT_LEAST_2_7)),cvsRoberts_ASAi_sps.ml) \
	   $(if $(and $(SUPERLUMT_ENABLED),\
	   	      $(AT_LEAST_2_7),\
		      $(OK_CVodeGetSensNumNonlinSolvIters)),cvsRoberts_FSA_sps.ml) \
	   $(if $(AT_LEAST_3_1_2),cvsRoberts_FSA_dns_Switch.ml) \
	   $(if $(AT_LEAST_3_2_0),cvsRoberts_FSA_dns_constraints.ml) \
	   $(if $(AT_LEAST_3_2_0),cvsRoberts_ASAi_dns_constraints.ml) \
	   $(if $(AT_LEAST_6_2_0),cvsParticle_dns.ml) \
	   $(if $(AT_LEAST_6_2_0),cvsPendulum_dns.ml)

FILES_TO_CLEAN = cvsRoberts_ASAi_dns_fwd_stats.csv \
		 cvsRoberts_ASAi_dns_bkw1_stats.csv \
		 cvsRoberts_ASAi_dns_bkw2_stats.csv \
		 cvsRoberts_FSA_dns_stats_-sensi_sim_t.csv \
		 cvsRoberts_FSA_dns_stats_-sensi_sim_f.csv \
		 cvsRoberts_FSA_dns_stats_-sensi_stg_t.csv \
		 cvsRoberts_FSA_dns_stats_-sensi_stg_f.csv \
		 cvsRoberts_FSA_dns_stats_-sensi_stg1_t.csv \
		 cvsRoberts_FSA_dns_stats_-sensi_stg1_f.csv \
		 cvsPendulum_dns_ref.txt \
		 cvsPendulum_dns_rtol_*.txt

include ../../examples.mk

test:
	echo $(SUNDIALS_VERSION)
	echo $(OK_CVodeGetSensNumNonlinSolvIters)
	echo "$(filter $(SUNDIALS_VERSION),4.0.0 4.0.1 4.0.2 4.1.0 5.0.0 5.1.0 5.2.0 5.3.0)"
# Tests with command-line arguments (make same choice as Sundials distrib)

$(eval $(call EXECUTION_RULE,cvsAdvDiff_FSA_non,$$< -sensi sim t))
$(eval $(call EXECUTION_RULE,cvsDiurnal_FSA_kry,$$< -sensi sim t))
$(eval $(call EXECUTION_RULE,cvsRoberts_FSA_dns,$$< -sensi sim t))
$(eval $(call EXECUTION_RULE,cvsRoberts_FSA_dns_constraints,$$< -sensi stg1 t))
$(eval $(call EXECUTION_RULE,cvsRoberts_FSA_klu,$$< -sensi sim t))
$(eval $(call EXECUTION_RULE,cvsRoberts_FSA_sps,$$< -sensi sim t))
