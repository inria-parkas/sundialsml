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
		     cvDirectDemo_ls.ml \
		     cvDiurnal_kry.ml \
		     cvDiurnal_kry_bp.ml \
		     cvKrylovDemo_ls.ml \
		     cvKrylovDemo_prec.ml \
		     $(if $(KLU_ENABLED),cvRoberts_klu.ml) \
		     $(if $(SUPERLUMT_ENABLED),cvRoberts_sps.ml) \

EXAMPLES_SINCE_3_0 =	cvDiurnal_kry_alt.ml \
 			cvRoberts_dns_uw_alt.ml
EXAMPLES_SINCE_3_1_2 =	cvDisc_dns.ml \
			cvRoberts_dns_negsol.ml

EXAMPLES_SINCE_3_2_0 = cvRoberts_dns_constraints.ml

EXAMPLES_SINCE_5_1_0 = cvRoberts_block_klu.ml

EXAMPLES_SINCE_5_3_0 = cvParticle_dns.ml \
		       cvPendulum_dns.ml

EXAMPLES = $(EXAMPLES_SINCE_2_5)				\
	   $(if $(AT_LEAST_3_0),$(EXAMPLES_SINCE_3_0))		\
	   $(if $(AT_LEAST_3_1_2),$(EXAMPLES_SINCE_3_1_2))	\
	   $(if $(AT_LEAST_3_2_0),$(EXAMPLES_SINCE_3_2_0))	\
	   $(if $(AT_LEAST_5_1_0),$(EXAMPLES_SINCE_5_1_0))	\
	   $(if $(AT_LEAST_5_3_0),$(EXAMPLES_SINCE_5_3_0))

LAPACK_EXAMPLES = cvAdvDiff_bndL.ml \
		  cvRoberts_dnsL.ml

FILES_TO_CLEAN = cvParticle_solution.txt cvParticle_error.txt \
	 cvPendulum_dns_ref.txt \
	 cvPendulum_dns_tol_1.00e-05.txt cvPendulum_dns_tol_1.00e-05_proj.txt \
	 cvPendulum_dns_tol_1.00e-06.txt cvPendulum_dns_tol_1.00e-06_proj.txt \
	 cvPendulum_dns_tol_1.00e-07.txt cvPendulum_dns_tol_1.00e-07_proj.txt \
	 cvPendulum_dns_tol_1.00e-08.txt cvPendulum_dns_tol_1.00e-08_proj.txt \
	 cvPendulum_dns_tol_1.00e-09.txt cvPendulum_dns_tol_1.00e-09_proj.txt \

include ../../examples.mk

$(eval $(call EXECUTION_RULE,cvRoberts_dns_uw_alt,$$<,cvRoberts_dns_uw))
$(eval $(call EXECUTION_RULE,cvDiurnal_kry_alt,$$<,cvDiurnal_kry))

