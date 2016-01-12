# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of IDAS.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless the
# latter are recompiled to link against IDAS.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=ida/serial

EXAMPLES = idaRoberts_dns.ml  \
	   idaHeat2D_bnd.ml   \
	   idaHeat2D_kry.ml   \
	   idaFoodWeb_bnd.ml  \
	   idaSlCrank_dns.ml  \
	   idaKrylovDemo_ls.ml \
	   idaHeat2D_klu.ml \
	   idaRoberts_klu.ml \
	   idaRoberts_dns_alt.ml

KLU_EXAMPLES = idaHeat2D_klu.ml \
	       idaRoberts_klu.ml
SUPERLUMT_EXAMPLES = idaRoberts_sps.ml

include ../../examples.mk

$(eval $(call EXECUTION_RULE,idaRoberts_dns_alt,$$<,idaRoberts_dns))
