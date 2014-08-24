include config

### Objects shared between sundials.cma and sundials_nosensi.cma.

# Common to CVODE, IDA, and KINSOL.
COBJ_COMMON = sundials_ml$(XO) dls_ml$(XO) nvector_ml$(XO) spils_ml$(XO)

COBJ_MAIN = $(COBJ_COMMON) kinsol_ml$(XO)

MLOBJ_MAIN = sundials.cmo dls.cmo spils.cmo nvector.cmo			\
	     nvector_custom.cmo nvector_array.cmo nvector_serial.cmo 	\
	     cvode_impl.cmo ida_impl.cmo kinsol_impl.cmo		\
	     cvode.cmo kinsol.cmo ida.cmo

### Objects specific to sundials.cma.
COBJ_SENS  = cvode_ml_s$(XO) ida_ml_s$(XO) cvodes_ml.o idas_ml.o
MLOBJ_SENS = cvodes.cmo idas.cmo

### Objects specific to sundials_nosensi.cma.
COBJ_NOSENSI = cvode_ml$(XO) ida_ml$(XO)
MLBJ_NOSENSI = 

### Objects specific to sundials_mpi.cma.
COBJ_MPI = nvector_parallel_ml.o kinsol_bbd_ml.o		\
	   cvode_bbd_ml.o cvodes_bbd_ml.o			\
	   ida_bbd_ml.o idas_bbd_ml.o
MLOBJ_MPI = nvector_parallel.cmo kinsol_bbd.cmo	\
	    cvode_bbd.cmo cvodes_bbd.cmo	\
	    ida_bbd.cmo idas_bbd.cmo

MPI_LIBLINK= -lsundials_nvecparallel

### Other sets of files.

# For `make clean'.  All object files, including ones that may not be
# built/updated under the current configuration.  Duplicates OK.
ALL_COBJ = $(COBJ_MAIN) $(COBJ_SENS) $(COBJ_NOSENSI) $(COBJ_MPI)
ALL_MLOBJ = $(MLOBJ_MAIN) $(MLOBJ_SENS) $(MLOBJ_NOSENSI) $(MLOBJ_MPI)
ALL_CMA = sundials.cma sundials_nosensi.cma sundials_mpi.cma

# Installed files.

INSTALL_CMA=sundials.cma sundials_nosensi.cma \
	    $(if $(MPI_ENABLED), sundials_mpi.cma)

STUBLIBS=$(foreach file,$(INSTALL_CMA:.cma=$(XS)), dllml$(file))

INSTALL_FILES=			\
    META			\
    $(MLOBJ_MAIN:.cmo=.cmi)	\
    $(MLOBJ_SENS:.cmo=.cmi)	\
    $(MLOBJ_NOSENSI:.cmo=.cmi)	\
    $(MLOBJ_MPI:.cmo=.cmi)	\
    $(INSTALL_CMA)		\
    $(INSTALL_CMA:.cma=.cmxa)	\
    $(INSTALL_CMA:.cma=$(XA))	\
    $(foreach file,$(INSTALL_CMA:.cma=$(XA)), libml$(file))

### Build rules.

.PHONY: all sundials install doc clean distclean

all: $(INSTALL_CMA) $(INSTALL_CMA:.cma=.cmxa)

sundials.cma sundials.cmxa: $(MLOBJ_MAIN) $(MLOBJ_SENS)			    \
			    $(MLOBJ_MAIN:.cmo=.cmx) $(MLOBJ_SENS:.cmo=.cmx) \
			    $(COBJ_MAIN) $(COBJ_SENS)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)	\
	    -o sundials -oc mlsundials $^	\
	    $(LIB_PATH)				\
	    $(OCAML_CVODES_LIBLINK)		\
	    $(OCAML_IDAS_LIBLINK)		\
	    $(OCAML_KINSOL_LIBLINK)

sundials_nosensi.cma sundials_nosensi.cmxa:				  \
			$(MLOBJ_MAIN) $(MLOBJ_NOSENSI)			  \
			$(MLOBJ_MAIN:.cmo=.cmx) $(MLOBJ_NOSENSI:.cmo=.cmx) \
			$(COBJ_MAIN) $(COBJ_NOSENSI)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)			\
	    -o sundials_nosensi -oc mlsundials_nosensi $^	\
	    $(LIB_PATH)						\
	    $(OCAML_CVODE_LIBLINK)				\
	    $(OCAML_IDA_LIBLINK)				\
	    $(OCAML_KINSOL_LIBLINK)

sundials_mpi.cma sundials_mpi.cmxa: $(MLOBJ_MPI) $(MLOBJ_MPI:.cmo=.cmx) \
				    $(COBJ_MPI)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)		\
	    -o sundials_mpi -oc mlsundials_mpi $^	\
	    $(LIB_PATH) $(MPI_LIBLINK)

# The CFLAGS settings for CVODE works for modules common to CVODE and IDA.
$(COBJ_COMMON): %.o: %.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

nvector_parallel_ml.o: nvector_parallel_ml.c
	$(MPICC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

# KINSOL-specific C files.
kinsol_ml.o: kinsol_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

kinsol_bbd_ml.o: kinsol_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

# CVODE[S]-specific C files.
cvode_ml.o: cvode_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvode_ml_s.o: cvode_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) \
	    -o $@ -c $<

cvodes_ml.o: cvodes_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) -o $@ -c $<

cvode_bbd_ml.o: cvode_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvodes_bbd_ml.o: cvodes_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) -o $@ -c $<

# IDA[S]-specific C files.
ida_ml.o: ida_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

ida_ml_s.o: ida_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) $(IDA_CFLAGS) \
	    -o $@ -c $<

idas_ml.o: idas_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) -o $@ -c $<

ida_bbd_ml.o: ida_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

idas_bbd_ml.o: idas_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) -o $@ -c $<

# Docs.
dochtml.cmo: INCLUDES += -I +ocamldoc
dochtml.cmo: OCAMLFLAGS += -pp "cpp $(CPPFLAGS) -DOCAML_3X=$(OCAML_3X)"

META: META.in
	@$(ECHO) "version = \"$(VERSION)\"" > $@
	@$(CAT) $< >> $@

doc: doc/html/index.html

doc/html/index.html: doc/html dochtml.cmo intro.doc 			\
		     $(filter-out %_impl.cmi, $(MLOBJ_MAIN:.cmo=.cmi))	\
		     $(MLOBJ_SENS:.cmo=.cmi) 				\
		     $(if $(MPI_ENABLED), $(MLOBJ_MPI:.cmo=.cmi))
	$(OCAMLDOC) -g dochtml.cmo $(INCLUDES) 			\
	    -cvode-doc-root "$(CVODE_DOC_ROOT)" 		\
	    -cvodes-doc-root "$(CVODES_DOC_ROOT)" 		\
	    -ida-doc-root "$(IDA_DOC_ROOT)" 			\
	    -idas-doc-root "$(IDAS_DOC_ROOT)" 			\
	    -kinsol-doc-root "$(KINSOL_DOC_ROOT)" 		\
	    -pp "$(DOCPP)"					\
	    -d ./doc/html/					\
	    -hide Cvode_impl,Ida_impl,Kinsol_impl		\
	    -t "Sundials"					\
	    -intro intro.doc					\
	    $(filter-out %_impl.mli, $(MLOBJ_MAIN:.cmo=.mli))	\
	    $(if $(MPI_ENABLED), $(MLOBJ_MPI:.cmo=.mli))	\
	    $(MLOBJ_SENS:.cmo=.mli)

doc/html:
	mkdir $@

### Install / Uninstall

install: $(INSTALL_CMA) $(INSTALL_CMA:.cma=.cmxa) doc META
	$(MKDIR) $(PKGDIR)
	$(CP) $(INSTALL_FILES) $(PKGDIR)
	$(CP) $(STUBLIBS) $(STUBDIR)
ifeq ($(INSTALL_DOCS), 1)
	$(MKDIR) $(DOCDIR)/html
	$(CP) doc/html/style.css doc/html/*.html $(DOCDIR)/html/
endif

uninstall:
	for f in $(STUBLIBS); do	 \
	    $(RM) $(STUBDIR)$$f || true; \
	done
	for f in $(INSTALL_FILES); do	 \
	    $(RM) $(PKGDIR)$$f || true;  \
	done
	-$(RMDIR) $(PKGDIR)
ifeq ($(INSTALL_DOCS), 1)
	-$(RM) $(DOCDIR)/html/style.css $(DOCDIR)/html/*.html
	-$(RMDIR) $(DOCDIR)/html
	-$(RMDIR) $(DOCDIR)
endif

ocamlfind: $(INSTALL_CMA) $(INSTALL_CMA:.cma=.cmxa) META
	ocamlfind install sundials $(INSTALL_FILES) $(STUBLIBS)

### Misc

depend: .depend
.depend:
	$(OCAMLDEP) \
	    -pp "cpp $(CPPFLAGS) -DOCAML_3X=$(OCAML_3X)" \
	    *.mli *.ml > .depend
	$(CC) -MM $(CFLAGS) *.c >> .depend

clean:
	-@($(MAKE) -C examples clean)
	-@$(RM) -f $(ALL_MLOBJ) $(ALL_MLOBJ:.cmo=.cmx) $(ALL_MLOBJ:.cmo=.o)
	-@$(RM) -f $(ALL_MLOBJ:.cmo=.cmi) $(ALL_MLOBJ:.cmo=.annot) $(ALL_COBJ)
	-@$(RM) -f $(ALL_CMA) $(ALL_CMA:.cma=.cmxa) $(ALL_CMA:.cma=.a)
	-@$(RM) -f $(foreach file,$(INSTALL_CMA:.cma=$(XA)),libml$(file))
	-@$(RM) -f $(foreach file,$(INSTALL_CMA:.cma=$(XS)),dllml$(file))
	-@$(RM) -f $(STUBLIBS)
	-@$(RM) -f dochtml.cmi dochtml.cmo

cleandoc:
	-@$(RM) -f doc/html/*.html doc/html/style.css

distclean: clean cleandoc
	-@($(MAKE) -C examples distclean)
	-@$(RM) -f META
	-@$(RM) -f config config.h

-include .depend
