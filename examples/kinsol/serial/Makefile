# Invoke `make' with USELIB=sundials to run the tests with the
# sundials.cma rather than sundials_no_sens.cma.  It should make no
# difference to the test results.
USELIB ?= sundials_no_sens
SRCROOT = ../../..
SUBDIR = kinsol/serial

EXAMPLES = kinFerTron_dns.ml			\
	   kinFoodWeb_kry.ml			\
	   kinFoodWeb_kry_custom.ml		\
	   kinKrylovDemo_ls.ml			\
	   kinLaplace_bnd.ml			\
	   kinLaplace_picard_bnd.ml		\
	   kinRoberts_fp.ml			\
	   kinRoboKin_dns.ml			\
	   kinRoboKin_dns_alt.ml

include ../../examples.mk

$(eval $(call EXECUTION_RULE,kinRoboKin_dns_alt,$$<,kinRoboKin_dns))
$(eval $(call EXECUTION_RULE,kinFoodWeb_kry_custom,$$<,kinFoodWeb_kry))
