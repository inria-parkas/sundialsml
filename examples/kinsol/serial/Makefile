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

kinRoboKin_dns_alt.sundials.out: kinRoboKin_dns.sundials.out
	cp $< $@
kinFoodWeb_kry_custom.sundials.out: kinFoodWeb_kry.sundials.out
	cp $< $@

kinRoboKin_dns_alt.sundials.time: kinRoboKin_dns.sundials.time
	cp $< $@
kinFoodWeb_kry_custom.sundials.time: kinFoodWeb_kry.sundials.time
	cp $< $@
