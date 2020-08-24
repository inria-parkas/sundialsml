# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=arkode/C_serial

EXAMPLES_SINCE_2_5 = ark_KrylovDemo_prec.ml \
		     ark_analytic.ml \
		     ark_brusselator.ml \
		     ark_brusselator_fp.ml \
		     ark_brusselator1D.ml \
		     ark_heat1D.ml \
		     ark_heat1D_adapt.ml \
		     ark_robertson.ml \
		     ark_robertson_root.ml \
		     $(if $(KLU_ENABLED),ark_brusselator1D_klu.ml) \
		     $(if $(SUPERLUMT_ENABLED),ark_brusselator1D_FEM_slu.ml)

EXAMPLES_SINCE_4_0_0 = ark_analytic_nonlin.ml

EXAMPLES = $(EXAMPLES_SINCE_2_5)				\
	   $(if $(AT_LEAST_4_0_0),$(EXAMPLES_SINCE_4_0_0))

FILES_TO_CLEAN = bruss_mesh.txt bruss_u.txt bruss_v.txt bruss_w.txt \
		 heat1D.txt heat_mesh.txt solution.txt

include ../../examples.mk

#TODO:
#$(eval $(call EXECUTION_RULE,cvRoberts_dns_uw_alt,$$<,cvRoberts_dns_uw))
