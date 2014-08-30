# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_nosensi
SRCROOT=../../..
SUBDIR=cvode/parallel

MPI_EXAMPLES = cvDiurnal_kry_p.ml \
	       cvDiurnal_kry_bbd_p.ml \
	       cvAdvDiff_non_p.ml

include ../../examples.mk

$(eval $(call EXECUTION_RULE,cvDiurnal_kry_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,cvDiurnal_kry_bbd_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,cvAdvDiff_non_p,$(MPIRUN) -np 6 $$<))
