SRCROOT = ../../..
SUBDIR = idas/serial

EXAMPLES = idasRoberts_FSA_dns.ml	\
	   idasAkzoNob_ASAi_dns.ml	\
	   idasAkzoNob_dns.ml		\
           idasHessian_ASA_FSA.ml	\
	   idasRoberts_ASAi_dns.ml	\
	   idasSlCrank_FSA_dns.ml	\
	   $(if $(and $(AT_LEAST_2_7),$(KLU_ENABLED)),idasRoberts_ASAi_klu.ml) \
	   $(if $(and $(AT_LEAST_2_7),$(KLU_ENABLED)),idasRoberts_FSA_klu.ml) \
	   $(if $(and $(AT_LEAST_2_7),$(SUPERLUMT_ENABLED)),idasRoberts_ASAi_sps.ml) \
	   $(if $(and $(AT_LEAST_2_7),$(SUPERLUMT_ENABLED)),idasRoberts_FSA_sps.ml)

include ../../examples.mk

# Tests with command-line arguments (make same choice as Sundials distrib)

$(eval $(call EXECUTION_RULE,idasRoberts_FSA_dns,$$< -sensi sim t))
$(eval $(call EXECUTION_RULE,idasRoberts_FSA_klu,$$< -sensi sim t))
$(eval $(call EXECUTION_RULE,idasRoberts_FSA_sps,$$< -sensi sim t))
