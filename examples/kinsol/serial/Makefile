# Invoke `make' with USELIB=sundials to run the tests with the
# sundials.cma rather than sundials_nosensi.cma.  It should make no
# difference to the test results.
USELIB ?= sundials_nosensi
SRCROOT = ../../..
SUBDIR = kinsol/serial

EXAMPLES = kinFerTron_dns.ml			\
	   kinFoodWeb_kry.ml			\
	   kinFoodWeb_kry_custom.ml		\
	   kinKrylovDemo_ls.ml			\
	   kinLaplace_bnd.ml			\
	   kinRoboKin_dns.ml			\
	   kinRoboKin_dns_alt.ml

include ../../examples.mk

$(eval $(call C_ALIAS,kinRoboKin_dns_alt,kinRoboKin_dns))
$(eval $(call C_ALIAS,kinFoodWeb_kry_custom,kinFoodWeb_kry))
