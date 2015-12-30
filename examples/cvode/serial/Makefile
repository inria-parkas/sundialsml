# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=cvode/serial

EXAMPLES = cvAdvDiff_bnd.ml \
	   cvRoberts_dns.ml \
	   cvRoberts_dns_uw.ml \
	   cvRoberts_dns_uw_alt.ml \
	   cvDirectDemo_ls.ml \
	   cvDiurnal_kry.ml \
	   cvDiurnal_kry_bp.ml \
	   cvKrylovDemo_ls.ml \
	   cvKrylovDemo_prec.ml

LAPACK_EXAMPLES = cvAdvDiff_bndL.ml \
		  cvRoberts_dnsL.ml

KLU_EXAMPLES = cvRoberts_klu.ml
SUPERLUMT_EXAMPLES = cvRoberts_sps.ml

include ../../examples.mk

$(eval $(call EXECUTION_RULE,cvRoberts_dns_uw_alt,$$<,cvRoberts_dns_uw))
