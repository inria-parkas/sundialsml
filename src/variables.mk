include ../config

### Common rules
# We can't use pattern rules here because they don't imply
# dependencies and ocamldep doesn't generate dependencies of the form
# foo.cm[iox]: foo.ml or foo.cmi: foo.mli
.SUFFIXES : .mli .ml .cmi .cmo .cmx

.ml.cmo:
	$(OCAMLC) $(OCAMLFLAGS) $(SUBDIRS:%=-I %) -c $(INCLUDES) $<

.mli.cmi:
	$(OCAMLC) $(OCAMLFLAGS) $(SUBDIRS:%=-I %) -c $(INCLUDES) $<

.ml.cmx:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) $(SUBDIRS:%=-I %) -c $(INCLUDES) $<

%.o: %.c
	$(CC) -I $(OCAML_INCLUDE) $(CFLAGS) $(CSUBDIRS) -o $@ -c $<

### Objects shared between sundials.cma and sundials_no_sens.cma.

# Common to CVODE, IDA, KINSOL, and ARKODE.
COBJ_COMMON = sundials/sundials_ml$(XO)	\
	      lsolvers/dls_ml$(XO)	\
	      $(SLS_ML_XO)		\
	      nvectors/nvector_ml$(XO)	\
	      lsolvers/spils_ml$(XO)	\
	      $(NVECPTHREADS_ML_XO)	\
	      $(NVECOPENMP_ML_XO)

COBJ_MAIN = $(COBJ_COMMON) kinsol/kinsol_ml$(XO) $(ARKODE_COBJ_MAIN)

MLOBJ_MAIN =	sundials/sundials_config.cmo	\
		sundials/sundials.cmo		\
		nvectors/nvector.cmo		\
		lsolvers/dls_impl.cmo		\
		lsolvers/dls.cmo		\
		lsolvers/sls_impl.cmo		\
		$(SLS_CMO)			\
		lsolvers/spils.cmo		\
		nvectors/nvector_custom.cmo	\
		nvectors/nvector_array.cmo	\
		nvectors/nvector_serial.cmo	\
		$(NVECPTHREADS_CMO)		\
		$(NVECOPENMP_CMO)		\
		cvode/cvode_impl.cmo		\
		ida/ida_impl.cmo		\
		kinsol/kinsol_impl.cmo		\
		cvode/cvode.cmo			\
		kinsol/kinsol.cmo		\
		ida/ida.cmo			\
		$(ARKODE_MLOBJ_MAIN)		\
		$(KLU_MLOBJ_MAIN)		\
		$(SUPERLUMT_MLOBJ_MAIN)

CMI_MAIN = $(filter-out sundials/sundials_config.cmi,$(filter-out %_impl.cmi,\
	    $(MLOBJ_MAIN:.cmo=.cmi)))

### Objects specific to sundials.cma.
COBJ_SENS  =	cvodes/cvode_ml_s$(XO)		\
		idas/ida_ml_s$(XO)		\
		cvodes/cvodes_ml.o		\
		idas/idas_ml.o			\
		$(KLU_COBJ_SENS)		\
		$(SUPERLUMT_COBJ_SENS)
MLOBJ_SENS =	cvodes/cvodes.cmo		\
		idas/idas.cmo			\
		$(KLU_MLOBJ_SENS)		\
		$(SUPERLUMT_MLOBJ_SENS)
CMI_SENS = $(MLOBJ_SENS:.cmo=.cmi)

### Objects specific to sundials_no_sens.cma.
COBJ_NO_SENS =	cvode/cvode_ml$(XO)		\
		ida/ida_ml$(XO)			\
		$(KLU_COBJ_NO_SENS)		\
		$(SUPERLUMT_COBJ_NO_SENS)
MLOBJ_NO_SENS =

### Objects specific to sundials_mpi.cma.
COBJ_MPI =	nvectors/nvector_parallel_ml$(XO)	\
		kinsol/kinsol_bbd_ml$(XO)		\
		$(ARKODE_COBJ_BBD)			\
		cvode/cvode_bbd_ml$(XO)			\
		cvodes/cvodes_bbd_ml$(XO)		\
		ida/ida_bbd_ml$(XO)			\
		idas/idas_bbd_ml$(XO)
MLOBJ_MPI =	nvectors/nvector_parallel.cmo	\
		kinsol/kinsol_bbd.cmo		\
		$(ARKODE_MLOBJ_BBD)		\
		cvode/cvode_bbd.cmo		\
		cvodes/cvodes_bbd.cmo		\
		ida/ida_bbd.cmo			\
		idas/idas_bbd.cmo
CMI_MPI = $(MLOBJ_MPI:.cmo=.cmi)

### Other sets of files.

# For `make clean'.  All object files, including ones that may not be
# built/updated under the current configuration.  Duplicates OK.
ALL_COBJ = $(COBJ_MAIN) $(COBJ_SENS) $(COBJ_NO_SENS) $(COBJ_MPI)
ALL_MLOBJ =doc/dochtml.cmo $(MLOBJ_MAIN) \
	   $(MLOBJ_SENS) $(MLOBJ_NO_SENS) $(MLOBJ_MPI)
ALL_CMA = sundials.cma sundials_no_sens.cma sundials_mpi.cma \
          sundials_docs.cma sundials_docs.cmxs

# Installed files.

INSTALL_CMA=sundials.cma sundials_no_sens.cma \
	    $(if $(MPI_ENABLED),sundials_mpi.cma)

INSTALL_CMI=$(CMI_MAIN) $(CMI_SENS) $(if $(MPI_ENABLED),$(CMI_MPI))

STUBLIBS=$(foreach file,$(INSTALL_CMA:.cma=$(XS)), dllml$(file))

# Don't include $(STUBLIBS) here; they go in a different directory.
INSTALL_FILES=							\
    META							\
    $(INSTALL_CMI)						\
    $(INSTALL_CMA)						\
    $(INSTALL_CMA:.cma=.cmxa)					\
    $(INSTALL_CMA:.cma=$(XA))					\
    $(foreach file,$(INSTALL_CMA:.cma=$(XA)), libml$(file))
