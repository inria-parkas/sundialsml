include variables.mk

default: all

# To compile with profiling (with gcc):
# ./configure CFLAGS=-pg OCAMLOPTFLAGS=-p ...

### Build rules.

.PHONY: all sundials install doc clean distclean

all: $(INSTALL_CMA) $(INSTALL_CMXA)

sundials.cma sundials.cmxa: $(MLOBJ_MAIN) $(MLOBJ_SENS)			    \
			    $(MLOBJ_MAIN:.cmo=.cmx) $(MLOBJ_SENS:.cmo=.cmx) \
			    $(COBJ_MAIN) $(COBJ_SENS)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)	\
	    $(if $(ENABLE_SHARED),,-custom)	\
	    -o sundials -oc mlsundials $^	\
	    $(LIB_PATH)				\
	    $(OCAML_CVODES_LIBLINK)		\
	    $(OCAML_ARKODE_LIBLINK)		\
	    $(OCAML_IDAS_LIBLINK)		\
	    $(OCAML_KINSOL_LIBLINK)		\
	    $(OCAML_ALL_LIBLINK)
sundials.cma: | sundials.cmxa # prevent simultaneous builds

sundials_no_sens.cma sundials_no_sens.cmxa:				  \
			$(MLOBJ_MAIN) $(MLOBJ_NO_SENS)			  \
			$(MLOBJ_MAIN:.cmo=.cmx) $(MLOBJ_NO_SENS:.cmo=.cmx) \
			$(COBJ_MAIN) $(COBJ_NO_SENS)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)			\
	    $(if $(ENABLE_SHARED),,-custom)			\
	    -o sundials_no_sens -oc mlsundials_no_sens $^	\
	    $(LIB_PATH)						\
	    $(OCAML_CVODE_LIBLINK)				\
	    $(OCAML_ARKODE_LIBLINK)				\
	    $(OCAML_IDA_LIBLINK)				\
	    $(OCAML_KINSOL_LIBLINK)				\
	    $(OCAML_ALL_LIBLINK)
sundials_no_sens.cma: | sundials_no_sens.cmxa # prevent simultaneous builds

sundials_mpi.cma sundials_mpi.cmxa: $(MLOBJ_MPI) $(MLOBJ_MPI:.cmo=.cmx) \
				    $(COBJ_MPI)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)		\
	    $(if $(ENABLE_SHARED),,-custom)		\
	    -o sundials_mpi -oc mlsundials_mpi $^	\
	    $(LIB_PATH) $(MPI_LIBLINK)
sundials_mpi.cma: | sundials_mpi.cmxa # prevent simultaneous builds

sundials_pthreads.cmxa sundials_pthreads.cma: $(MLOBJ_PTHREADS)		  \
					      $(MLOBJ_PTHREADS:.cmo=.cmx) \
					      $(COBJ_PTHREADS)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)		\
	    $(if $(ENABLE_SHARED),,-custom)		\
	    -o sundials_pthreads -oc mlsundials_pthreads $^	\
	    $(LIB_PATH) $(PTHREADS_LIBLINK)
sundials_pthreads.cma: | sundials_pthreads.cmxa # prevent simultaneous builds

sundials_openmp.cma sundials_openmp.cmxa: $(MLOBJ_OPENMP)		\
					  $(MLOBJ_OPENMP:.cmo=.cmx)	\
					  $(COBJ_OPENMP)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)		\
	    $(if $(ENABLE_SHARED),,-custom)		\
	    -o sundials_openmp -oc mlsundials_openmp $^	\
	    $(LIB_PATH) $(OPENMP_LIBLINK)
sundials_openmp.cma: | sundials_openmp.cmxa # prevent simultaneous builds

$(CMA_TOP_ALL): %.cma:
	$(OCAMLC) -a -o $@ $^

$(CMA_TOP_ALL:.cma=.cmxa): %.cmxa:
	$(OCAMLOPT) -a -o $@ $^

sundials_top.cma: $(MLOBJ_TOP)
sundials_top_mpi.cma: $(MLOBJ_TOP_MPI)
sundials_top_openmp.cma: $(MLOBJ_TOP_OPENMP)
sundials_top_pthreads.cma: $(MLOBJ_TOP_PTHREADS)
sundials_top_findlib.cma: $(MLOBJ_TOP_FINDLIB)

sundials_top.cmxa: $(MLOBJ_TOP:.cmo=.cmx)
sundials_top_mpi.cmxa: $(MLOBJ_TOP_MPI:.cmo=.cmx)
sundials_top_openmp.cmxa: $(MLOBJ_TOP_OPENMP:.cmo=.cmx)
sundials_top_pthreads.cmxa: $(MLOBJ_TOP_PTHREADS:.cmo=.cmx)
sundials_top_findlib.cmxa: $(MLOBJ_TOP_FINDLIB:.cmo=.cmx)

$(MLOBJ_TOP_ALL) $(CMI_TOP) $(MLOBJ_TOP_ALL:.cmo=.cmx):	\
	INCLUDES += $(TOP_INCLUDES)

$(MLOBJ_TOP_FINDLIB:.cmo=.ml): %.ml: ../config
	echo "(* Automatically generated file - don't edit!  See the Makefile.  *)" > $@
	ocamlfind printconf destdir | sed -e 's/^/Topdirs.dir_directory \"/' \
	    -e 's/$$/\/sundialsml";;/' >> $@

$(MLOBJ_MPI) $(CMI_MPI) $(MLOBJ_MPI:.cmo=.cmx)	\
	     doc/html/index.html :		\
    INCLUDES += $(MPI_INCLUDES)

# Avoid cyclical dependencies for modules aliased in sundials/sundials.mli
# see: https://caml.inria.fr/pub/docs/manual-ocaml/depend.html#sec334
# (14.2: A typical Makefile)
# see: https://caml.inria.fr/pub/docs/manual-ocaml/extn.html
# (8.12: Type-level module aliases)
sundials/sundials.cmi: sundials/sundials.mli
	$(OCAMLC) $(OCAMLFLAGS) -no-alias-deps -w -49 $(SUBDIRS:%=-I %) -c $(INCLUDES) $<

sundials/sundials.cmo: sundials/sundials.ml
	$(OCAMLC) $(OCAMLFLAGS) -no-alias-deps -w -49 $(SUBDIRS:%=-I %) -c $(INCLUDES) $<

sundials/sundials.cmx: sundials/sundials.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -no-alias-deps -w -49 $(SUBDIRS:%=-I %) -c $(INCLUDES) $<

$(filter-out sundials/sundials.cmi,$(MLOBJ_AFTER_SUNDIALS:.cmo=.cmi)) : \
    sundials/sundials.cmi

# The CFLAGS settings for CVODE works for modules common to CVODE and IDA.
$(COBJ_COMMON): %.o: %.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

nvectors/nvector_parallel_ml.o: nvectors/nvector_parallel_ml.c \
		sundials/sundials_ml.h nvectors/nvector_ml.h \
		nvectors/nvector_parallel_ml.h
	$(MPICC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

nvectors/nvector_openmp_ml.o: nvectors/nvector_openmp_ml.c \
		sundials/sundials_ml.h nvectors/nvector_ml.h \
		nvectors/nvector_openmp_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

nvectors/nvector_pthreads_ml.o: nvectors/nvector_pthreads_ml.c \
		sundials/sundials_ml.h nvectors/nvector_ml.h \
		nvectors/nvector_pthreads_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

# KINSOL-specific C files.
kinsol/kinsol_ml.o: kinsol/kinsol_ml.c \
		sundials/sundials_ml.h lsolvers/sundials_matrix_ml.h \
		lsolvers/sundials_linearsolver_ml.h kinsol/kinsol_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

kinsol/kinsol_bbd_ml.o: kinsol/kinsol_bbd_ml.c kinsol/kinsol_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

kinsol/kinsol_klu_ml.o: kinsol/kinsol_klu_ml.c \
		sundials/sundials_ml.h kinsol/kinsol_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

kinsol/kinsol_superlumt_ml.o: kinsol/kinsol_superlumt_ml.c \
		sundials/sundials_ml.h kinsol/kinsol_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

# ARKODE-specific C files.
arkode/arkode_ml.o: arkode/arkode_ml.c \
		lsolvers/sundials_matrix_ml.h \
		lsolvers/sundials_linearsolver_ml.h \
		sundials/sundials_ml.h arkode/arkode_ml.h nvectors/nvector_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) -o $@ -c $<

arkode/arkode_bbd_ml.o: arkode/arkode_bbd_ml.c \
		sundials/sundials_ml.h arkode/arkode_ml.h nvectors/nvector_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) -o $@ -c $<

arkode/arkode_klu_ml.o: arkode/arkode_klu_ml.c \
		sundials/sundials_ml.h arkode/arkode_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) -o $@ -c $<

arkode/arkode_superlumt_ml.o: arkode/arkode_superlumt_ml.c \
		sundials/sundials_ml.h arkode/arkode_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) -o $@ -c $<

# CVODE[S]-specific C files.
cvode/cvode_ml.o: cvode/cvode_ml.c \
		lsolvers/sundials_matrix_ml.h \
		lsolvers/sundials_linearsolver_ml.h \
		sundials/sundials_ml.h cvode/cvode_ml.h nvectors/nvector_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvodes/cvode_ml_s.o: cvode/cvode_ml.c \
		lsolvers/sundials_matrix_ml.h \
		lsolvers/sundials_linearsolver_ml.h \
		sundials/sundials_ml.h cvode/cvode_ml.h nvectors/nvector_ml.h
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(CVODE_CFLAGS) -o $@ -c $<

cvodes/cvodes_ml.o: cvodes/cvodes_ml.c \
		lsolvers/sundials_matrix_ml.h \
		lsolvers/sundials_linearsolver_ml.h \
		sundials/sundials_ml.h cvode/cvode_ml.h nvectors/nvector_ml.h \
		cvodes/cvodes_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) -Icvode -o $@ -c $<

cvode/cvode_bbd_ml.o: cvode/cvode_bbd_ml.c \
		sundials/sundials_ml.h cvode/cvode_ml.h nvectors/nvector_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvodes/cvodes_bbd_ml.o: cvodes/cvodes_bbd_ml.c \
		sundials/sundials_ml.h cvode/cvode_ml.h nvectors/nvector_ml.h \
		cvodes/cvodes_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) -Icvode -o $@ -c $<

cvode/cvode_klu_ml.o: cvode/cvode_klu_ml.c \
		sundials/sundials_ml.h cvode/cvode_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvodes/cvode_klu_ml_s.o: cvode/cvode_klu_ml.c \
		sundials/sundials_ml.h cvode/cvode_ml.h \
		lsolvers/sundials_matrix_ml.h cvodes/cvodes_ml.h
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(CVODE_CFLAGS) -o $@ -c $<

cvodes/cvodes_klu_ml.o: cvodes/cvodes_klu_ml.c \
		sundials/sundials_ml.h cvode/cvode_ml.h \
		lsolvers/sundials_matrix_ml.h cvodes/cvodes_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) -Icvode -o $@ -c $<

cvode/cvode_superlumt_ml.o: cvode/cvode_superlumt_ml.c \
		sundials/sundials_ml.h cvode/cvode_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvodes/cvode_superlumt_ml_s.o: cvode/cvode_superlumt_ml.c \
		sundials/sundials_ml.h cvode/cvode_ml.h \
		lsolvers/sundials_matrix_ml.h cvodes/cvodes_ml.h
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(CVODE_CFLAGS) -o $@ -c $<

cvodes/cvodes_superlumt_ml.o: cvodes/cvodes_superlumt_ml.c \
		sundials/sundials_ml.h cvode/cvode_ml.h \
		lsolvers/sundials_matrix_ml.h cvodes/cvodes_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) -Icvode -o $@ -c $<

# IDA[S]-specific C files.
ida/ida_ml.o: ida/ida_ml.c ida/ida_ml.h \
		sundials/sundials_ml.h lsolvers/sundials_linearsolver_ml.h \
		nvectors/nvector_ml.h lsolvers/sundials_matrix_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

idas/ida_ml_s.o: ida/ida_ml.c ida/ida_ml.h \
		sundials/sundials_ml.h lsolvers/sundials_linearsolver_ml.h \
		nvectors/nvector_ml.h lsolvers/sundials_matrix_ml.h
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(IDA_CFLAGS) -o $@ -c $<

idas/idas_ml.o: idas/idas_ml.c ida/ida_ml.h idas/idas_ml.h \
		lsolvers/sundials_linearsolver_ml.h \
		sundials/sundials_ml.h nvectors/nvector_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) -Iida -o $@ -c $<

ida/ida_bbd_ml.o: ida/ida_bbd_ml.c \
		sundials/sundials_ml.h ida/ida_ml.h nvectors/nvector_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

idas/idas_bbd_ml.o: idas/idas_bbd_ml.c lsolvers/sundials_matrix_ml.h \
		lsolvers/sundials_linearsolver_ml.h \
		sundials/sundials_ml.h ida/ida_ml.h idas/idas_ml.h \
		nvectors/nvector_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) -Iida -o $@ -c $<

ida/ida_klu_ml.o: ida/ida_klu_ml.c \
		sundials/sundials_ml.h ida/ida_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

idas/ida_klu_ml_s.o: ida/ida_klu_ml.c \
		sundials/sundials_ml.h ida/ida_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(IDA_CFLAGS) -o $@ -c $<

idas/idas_klu_ml.o: idas/idas_klu_ml.c \
		sundials/sundials_ml.h ida/ida_ml.h \
		lsolvers/sundials_matrix_ml.h idas/idas_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) -Iida -o $@ -c $<

ida/ida_superlumt_ml.o: ida/ida_superlumt_ml.c \
		sundials/sundials_ml.h ida/ida_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

idas/ida_superlumt_ml_s.o: ida/ida_superlumt_ml.c \
		sundials/sundials_ml.h ida/ida_ml.h \
		lsolvers/sundials_matrix_ml.h
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(IDA_CFLAGS) -o $@ -c $<

idas/idas_superlumt_ml.o: idas/idas_superlumt_ml.c \
		sundials/sundials_ml.h ida/ida_ml.h \
		lsolvers/sundials_matrix_ml.h idas/idas_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) -o $@ -c $<

# Other files

META: META.in ../config
	$(CPP) -DVERSION=\"$(VERSION)\"				\
	    $(if $(MPI_ENABLED),-DMPI_ENABLED)			\
	    $(if $(PTHREADS_ENABLED),-DPTHREADS_ENABLED)	\
	    $(if $(OPENMP_ENABLED),-DOPENMP_ENABLED)		\
	    $(if $(TOP_ENABLED),-DTOP_ENABLED)			\
	    $<							\
	    | grep -v '^#' > $@

### Install / Uninstall

install: install-sys $(if $(INSTALL_DOCS),install-doc)

# Install to OCaml's system directory -- /usr/lib/ocaml on Debian derivatives.
install-sys: META $(INSTALL_CMA) $(INSTALL_CMXA)
	[ -d $(PKGDIR) ] || $(MKDIR) $(PKGDIR)
	$(CP) $(INSTALL_FILES) $(PKGDIR)
	$(if $(ENABLE_SHARED),[ -d $(STUBDIR) ] || $(MKDIR) $(STUBDIR))
	$(if $(ENABLE_SHARED),$(CP) $(STUBLIBS) $(STUBDIR))

install-ocamlfind: install-findlib
install-findlib: META $(INSTALL_CMA) $(INSTALL_CMXA)		\
		 $(INSTALL_CMA_FINDLIB)
	@ocamlfind install sundialsml $(INSTALL_FILES)	\
	    $(INSTALL_FILES_FINDLIB)			\
	    -optional $(STUBLIBS)

install-doc:
	$(MAKE) -C ../doc install

uninstall-doc:
	$(MAKE) -C ../doc uninstall

uninstall: uninstall-sys uninstall-doc

uninstall-sys:
	-$(RM) $(foreach f,$(STUBLIBS),$(STUBDIR)$f)
	-$(RM) $(foreach f,$(notdir $(INSTALL_FILES)),$(PKGDIR)$f)
	-$(RMDIR) $(PKGDIR)

uninstall-ocamlfind: uninstall-findlib
uninstall-findlib:
	@ocamlfind remove sundialsml

### Misc

depend: .depend
.depend:
	$(OCAMLDEP) $(foreach d,$(SUBDIRS),-I $d) \
	    -map sundials/sundials.mli \
	    $(filter-out sundials/sundials.ml,$(ALL_MLOBJ:.cmo=.ml)) \
	    $(filter-out sundials/sundials.mli,$(ALL_MLOBJ:.cmo=.mli)) \
	    > .depend
	$(OCAMLDEP) $(foreach d,$(SUBDIRS),-I $d) \
	    -as-map sundials/sundials.mli sundials/sundials.ml >> .depend
	#$(CC) -MM $(CPPFLAGS) */*.c >> .depend

clean:
	-@$(RM) -f $(ALL_MLOBJ) $(ALL_MLOBJ:.cmo=.cmx) $(ALL_MLOBJ:.cmo=.o)
	-@$(RM) -f $(ALL_MLOBJ:.cmo=.cmi) $(ALL_MLOBJ:.cmo=.annot) $(ALL_COBJ)
	-@$(RM) -f $(ALL_MLOBJ:.cmo=.cmt) $(ALL_MLOBJ:.cmo=.cmti)
	-@$(RM) -f $(ALL_CMA) $(ALL_CMA:.cma=.cmxa) $(ALL_CMA:.cma=.a)
	-@$(RM) -f $(foreach file,$(ALL_CMA:.cma=$(XA)),libml$(file))
	-@$(RM) -f $(foreach file,$(ALL_CMA:.cma=$(XS)),dllml$(file))
	-@$(RM) -f $(STUBLIBS)

distclean: clean
	-@$(RM) -f META
	-@$(RM) -f sundials/sundials_Index.mli sundials/sundials_Index.ml
	-@$(RM) -f config.h sundials/sundials_configuration.ml

-include .depend
