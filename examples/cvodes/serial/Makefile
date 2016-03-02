SRCROOT=../../..
SUBDIR=cvodes/serial

EXAMPLES = cvsAdvDiff_ASAi_bnd.ml \
	   cvsAdvDiff_FSA_non.ml  \
	   cvsDiurnal_FSA_kry.ml  \
	   cvsFoodWeb_ASAi_kry.ml \
	   cvsFoodWeb_ASAp_kry.ml \
	   cvsHessian_ASA_FSA.ml  \
	   cvsRoberts_ASAi_dns.ml \
	   cvsRoberts_FSA_dns.ml  \
	   $(if $(KLU_ENABLED),$(KLU_EXAMPLES)) \
	   $(if $(SUPERLUMT_ENABLED),$(SUPERLUMT_EXAMPLES))

include ../../examples.mk

# Tests with command-line arguments (make same choice as Sundials distrib)

$(eval $(call EXECUTION_RULE,cvsAdvDiff_FSA_non,$$< -sensi sim t))
$(eval $(call EXECUTION_RULE,cvsDiurnal_FSA_kry,$$< -sensi sim t))
$(eval $(call EXECUTION_RULE,cvsRoberts_FSA_dns,$$< -sensi sim t))
