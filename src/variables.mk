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
	      lsolvers/spils_ml$(XO)

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
		cvode/cvode_impl.cmo		\
		ida/ida_impl.cmo		\
		kinsol/kinsol_impl.cmo		\
		cvode/cvode.cmo			\
		kinsol/kinsol.cmo		\
		ida/ida.cmo			\
		$(ARKODE_MLOBJ_MAIN)

CMI_MAIN = $(filter-out sundials/sundials_config.cmi,$(filter-out %_impl.cmi,\
	    $(MLOBJ_MAIN:.cmo=.cmi)))

### Objects specific to sundials.cma.
COBJ_SENS  =	cvodes/cvode_ml_s$(XO)		\
		idas/ida_ml_s$(XO)		\
		cvodes/cvodes_ml$(XO)		\
		idas/idas_ml$(XO)		\
		cvodes/cvode_klu_ml_s${XO}	\
		cvodes/cvode_superlumt_ml_s${XO}\
		cvodes/cvodes_klu_ml${XO}	\
		cvodes/cvodes_superlumt_ml${XO}	\
		arkode/arkode_klu_ml$(XO)	\
		arkode/arkode_superlumt_ml$(XO)	\
		kinsol/kinsol_klu_ml${XO}	\
		kinsol/kinsol_superlumt_ml${XO}	\
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
		cvode/cvode_klu_ml${XO}		\
		cvode/cvode_superlumt_ml${XO}	\
		ida/ida_klu_ml${XO}		\
		ida/ida_superlumt_ml${XO}	\
		cvodes/cvodes_klu_ml${XO}	\
		cvodes/cvodes_superlumt_ml${XO}	\
		arkode/arkode_klu_ml$(XO)	\
		arkode/arkode_superlumt_ml$(XO)	\
		kinsol/kinsol_klu_ml${XO}	\
		kinsol/kinsol_superlumt_ml${XO}
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

### Objects specific to sundials_pthreads.cma.
COBJ_PTHREADS =	nvectors/nvector_pthreads_ml$(XO)
MLOBJ_PTHREADS=	nvectors/nvector_pthreads.cmo
CMI_PTHREADS =	$(MLOBJ_PTHREADS:.cmo=.cmi)

### Objects specific to sundials_pthreads.cma.
COBJ_OPENMP =	nvectors/nvector_openmp_ml$(XO)
MLOBJ_OPENMP =	nvectors/nvector_openmp.cmo
CMI_OPENMP =	$(MLOBJ_OPENMP:.cmo=.cmi)

### Other sets of files.

# For `make clean'.  All object files, including ones that may not be
# built/updated under the current configuration.  Duplicates OK.
ALL_COBJ = $(COBJ_MAIN) $(COBJ_SENS) $(COBJ_NO_SENS) $(COBJ_MPI) \
	   $(COBJ_OPENMP) $(COBJ_PTHREADS)
ALL_MLOBJ =doc/dochtml.cmo $(MLOBJ_MAIN)		\
	   $(MLOBJ_SENS) $(MLOBJ_NO_SENS) $(MLOBJ_MPI)	\
	   $(MLOBJ_OPENMP) $(MLOBJ_PTHREADS)
ALL_CMA = sundials.cma sundials_no_sens.cma sundials_mpi.cma	\
	  sundials_openmp.cma sundials_pthreads.cma		\
	  sundials_docs.cma sundials_docs.cmxs

# Installed files.

INSTALL_CMA=sundials.cma sundials_no_sens.cma			\
	    $(if $(MPI_ENABLED),sundials_mpi.cma)		\
	    $(if $(OPENMP_ENABLED),sundials_openmp.cma)		\
	    $(if $(PTHREADS_ENABLED),sundials_pthreads.cma)

INSTALL_CMI=$(CMI_MAIN) $(CMI_SENS)			\
	    $(if $(MPI_ENABLED),$(CMI_MPI))		\
	    $(if $(PTHREADS_ENABLED),$(CMI_PTHREADS))	\
	    $(if $(OPENMP_ENABLED),$(CMI_OPENMP))

STUBLIBS=$(foreach file,$(INSTALL_CMA:.cma=$(XS)), dllml$(file))

# Don't include $(STUBLIBS) here; they go in a different directory.
INSTALL_FILES=							\
    META							\
    $(INSTALL_CMI)						\
    $(INSTALL_CMA)						\
    $(INSTALL_CMA:.cma=.cmxa)					\
    $(INSTALL_CMA:.cma=$(XA))					\
    $(foreach file,$(INSTALL_CMA:.cma=$(XA)), libml$(file))
