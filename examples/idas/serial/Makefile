SRCROOT = ../../..
SUBDIR = idas/serial

EXAMPLES = idasRoberts_FSA_dns.ml	\
	   idasAkzoNob_ASAi_dns.ml	\
	   idasAkzoNob_dns.ml		\
           idasHessian_ASA_FSA.ml	\
	   idasRoberts_ASAi_dns.ml	\
	   idasSlCrank_FSA_dns.ml

include ../../examples.mk

# Tests with command-line arguments (make same choice as Sundials distrib)

$(eval $(call EXECUTION_RULE,idasRoberts_FSA_dns,$$< -sensi sim t))
