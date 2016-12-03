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
	   idaRoberts_dns_alt.ml \
	   $(if $(KLU_ENABLED),idaHeat2D_klu.ml) \
	   $(if $(KLU_ENABLED),idaRoberts_klu.ml) \
	   $(if $(SUPERLUMT_ENABLED),idaRoberts_sps.ml) \
	   $(if $(AT_LEAST_2_7),idaFoodWeb_kry.ml)

include ../../examples.mk

$(eval $(call EXECUTION_RULE,idaRoberts_dns_alt,$$<,idaRoberts_dns))
