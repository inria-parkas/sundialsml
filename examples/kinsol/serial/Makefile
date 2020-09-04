# Invoke `make' with USELIB=sundials to run the tests with the
# sundials.cma rather than sundials_no_sens.cma.  It should make no
# difference to the test results.
USELIB ?= sundials_no_sens
SRCROOT = ../../..
SUBDIR = kinsol/serial

EXAMPLES_SINCE_2_5 = kinFerTron_dns.ml		\
		     kinFoodWeb_kry.ml		\
		     kinFoodWeb_kry_custom.ml	\
		     kinKrylovDemo_ls.ml	\
		     kinLaplace_bnd.ml		\
		     kinRoboKin_dns.ml		\
		     $(if $(and $(AT_LEAST_2_7),$(KLU_ENABLED)),kinFerTron_klu.ml)

EXAMPLES_SINCE_2_6 = kinRoberts_fp.ml		\
		     kinLaplace_picard_bnd.ml   \
		     $(if $(SUPERLUMT_ENABLED),kinRoboKin_slu.ml)

EXAMPLES_SINCE_3_0 = kinRoboKin_dns_alt.ml

EXAMPLES = $(EXAMPLES_SINCE_2_5) $(EXAMPLES_SINCE_2_6) $(EXAMPLES_SINCE_3_0)


ENABLED_EXAMPLES = $(EXAMPLES_SINCE_2_5)			\
		   $(if $(AT_LEAST_2_6),$(EXAMPLES_SINCE_2_6))	\
	   	   $(if $(AT_LEAST_3_0),$(EXAMPLES_SINCE_3_0))

include ../../examples.mk

$(eval $(call EXECUTION_RULE,kinRoboKin_dns_alt,$$<,kinRoboKin_dns))
$(eval $(call EXECUTION_RULE,kinFoodWeb_kry_custom,$$<,kinFoodWeb_kry))
