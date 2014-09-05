# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of IDAS.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless the
# latter are recompiled to link against IDAS.
USELIB ?= sundials_nosensi
SRCROOT=../../..
SUBDIR=ida/serial

EXAMPLES = idaRoberts_dns.ml  \
	   idaHeat2D_bnd.ml   \
	   idaHeat2D_kry.ml   \
	   idaFoodWeb_bnd.ml  \
	   idaSlCrank_dns.ml  \
	   idaKrylovDemo_ls.ml \
	   idaRoberts_dns_alt.ml

include ../../examples.mk

idaRoberts_dns_alt.sundials.out: idaRoberts_dns.sundials.out
	cp $< $@

idaRoberts_dns_alt.sundials.time: idaRoberts_dns.sundials.time
	cp $< $@
