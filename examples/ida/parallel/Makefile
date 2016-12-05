# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of IDAS.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless the
# latter are recompiled to link against IDAS.
USELIB ?= sundials_no_sens
SRCROOT=../../..
SUBDIR=ida/parallel
MPI_EXAMPLES = idaFoodWeb_kry_p.ml	\
	       idaFoodWeb_kry_bbd_p.ml	\
	       $(if $(AT_LEAST_2_7),idaHeat2D_kry_p.ml)	\
	       idaHeat2D_kry_bbd_p.ml

include ../../examples.mk

# Tests with MPI.

$(eval $(call EXECUTION_RULE,idaFoodWeb_kry_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,idaFoodWeb_kry_bbd_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,idaHeat2D_kry_bbd_p,$(MPIRUN) -np 4 $$<))
$(eval $(call EXECUTION_RULE,idaHeat2D_kry_p,$(MPIRUN) -np 4 $$<))
