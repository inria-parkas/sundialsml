# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=cvode/serial

EXAMPLES_SINCE_2_5 = cvAdvDiff_bnd.ml \
		     cvRoberts_dns.ml \
		     cvRoberts_dns_uw.ml \
		     cvRoberts_dns_uw_alt.ml \
		     cvDirectDemo_ls.ml \
		     cvDiurnal_kry.ml \
		     cvDiurnal_kry_bp.ml \
		     cvKrylovDemo_ls.ml \
		     cvKrylovDemo_prec.ml \
		     $(if $(KLU_ENABLED),cvRoberts_klu.ml) \
		     $(if $(SUPERLUMT_ENABLED),cvRoberts_sps.ml) \

EXAMPLES_SINCE_3_0 = cvDiurnal_kry_alt.ml

EXAMPLES_SINCE_3_1_2 = cvDisc_dns.ml

EXAMPLES = $(EXAMPLES_SINCE_2_5)				      \
	   $(if $(AT_LEAST_3_0),$(EXAMPLES_SINCE_3_0))		      \
	   $(if $(AT_LEAST_3_1_2),$(EXAMPLES_SINCE_3_1_2))

LAPACK_EXAMPLES = cvAdvDiff_bndL.ml \
		  cvRoberts_dnsL.ml

include ../../examples.mk

$(eval $(call EXECUTION_RULE,cvRoberts_dns_uw_alt,$$<,cvRoberts_dns_uw))
$(eval $(call EXECUTION_RULE,cvDiurnal_kry_alt,$$<,cvDiurnal_kry))

