include config

MLOBJ_MAIN = sundials.cmo dls.cmo \
	     nvector_serial.cmo nvector_custom.cmo nvector_array.cmo \
	     $(MPI_MODULES) \
	     spils.cmo cvode.cmo kinsol.cmo \
	     ida.cmo

MLOBJ_SENS = cvodes.cmo idas.cmo

MLOBJ_LOCAL = cvode_impl.cmo ida_impl.cmo kinsol_impl.cmo

MLOBJ_WOS = sundials.cmo $(MLOBJ_LOCAL) $(filter-out sundials.cmo,$(MLOBJ_MAIN))
MLOBJ = $(MLOBJ_WOS) $(MLOBJ_SENS)

COMMON_COBJ= sundials_ml$(XO) dls_ml$(XO) nvector_ml$(XO) spils_ml$(XO) \
	     $(MPI_COBJ)

COBJ_WOS = $(COMMON_COBJ) cvode_ml$(XO) ida_ml$(XO) kinsol_ml$(XO)
COBJ = $(COMMON_COBJ) cvode_ml_s$(XO) cvodes_ml$(XO) \
	ida_ml_s$(XO) idas_ml$(XO) kinsol_ml$(XO)

ALL_COBJ= $(COBJ) cvode_ml$(XO) ida_ml$(XO)

INSTALL_FILES=			\
    META			\
    $(MLOBJ_MAIN:.cmo=.cmi)	\
    $(MLOBJ_SENS:.cmo=.cmi)	\
    libmlsundials$(XA)		\
    sundials$(XA)		\
    sundials.cma		\
    sundials.cmxa		\
    libmlsundials_wos$(XA)	\
    sundials_wos$(XA)		\
    sundials_wos.cma		\
    sundials_wos.cmxa

STUBLIBS=dllmlsundials$(XS)

CFLAGS+=-fPIC

# ##

.PHONY: all sundials install doc

all: sundials.cma sundials.cmxa sundials_wos.cma sundials_wos.cmxa

# TODO: fix this:
sundials.cma sundials.cmxa: $(MLOBJ) $(MLOBJ:.cmo=.cmx) $(COBJ)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)	\
	    -o sundials -oc mlsundials $^	\
	    $(OCAML_CVODES_LIBLINK)		\
	    $(OCAML_IDAS_LIBLINK)		\
	    $(OCAML_KINSOL_LIBLINK)		\
	    $(NVECTOR_LIB)

# wos = without sensitivity
# TODO: fix this:
sundials_wos.cma sundials_wos.cmxa: $(MLOBJ_WOS) $(MLOBJ_WOS:.cmo=.cmx) \
				    $(COBJ_WOS)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)	  \
	    -o sundials_wos -oc mlsundials_wos $^ \
	    $(OCAML_CVODE_LIBLINK)		  \
	    $(OCAML_IDA_LIBLINK)		  \
	    $(OCAML_KINSOL_LIBLINK)		  \
	    $(NVECTOR_LIB)

# There are three sets of flags:
#   - one for CVODE-specific files
#   - one for IDA-specific files
#   - one for files common to CVODE and IDA

# The CFLAGS settings for CVODE works for modules common to CVODE and IDA.
sundials_ml.o: sundials_ml.c sundials_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

dls_ml.o: dls_ml.c dls_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

nvector_ml.o: nvector_ml.c nvector_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

nvector_parallel_ml.o: nvector_parallel_ml.c nvector_parallel_ml.h nvector_ml.h
	$(MPICC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvode_ml.o: cvode_ml.c dls_ml.h spils_ml.h cvode_ml.h sundials_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvode_ml_s.o: cvode_ml.c dls_ml.h spils_ml.h cvode_ml.h sundials_ml.h
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) \
	    -o $@ -c $<

cvodes_ml.o: cvodes_ml.c dls_ml.h spils_ml.h \
	     cvode_ml.h cvodes_ml.h sundials_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) -o $@ -c $<

ida_ml.o: ida_ml.c dls_ml.h spils_ml.h ida_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

ida_ml_s.o: ida_ml.c dls_ml.h spils_ml.h ida_ml.h
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) $(IDA_CFLAGS) \
	    -o $@ -c $<

idas_ml.o: idas_ml.c dls_ml.h spils_ml.h \
	     ida_ml.h idas_ml.h sundials_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) -o $@ -c $<

kinsol_ml.o: kinsol_ml.c dls_ml.h spils_ml.h kinsol_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

spils_ml.o: spils_ml.c sundials_ml.h spils_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

dochtml.cmo: INCLUDES += -I +ocamldoc
dochtml.cmo: OCAMLFLAGS += -pp "cpp $(CPPFLAGS) -DOCAML_3X=$(OCAML_3X)"

META: META.in
	@$(ECHO) "version = \"$(VERSION)\"" > $@
	@$(CAT) $< >> $@

doc: doc/html/index.html

doc/html/index.html: doc/html dochtml.cmo intro.doc \
		     $(MLOBJ_MAIN:.cmo=.mli) $(MLOBJ_MAIN:.cmo=.cmi)  \
		     $(MLOBJ_SENS:.cmo=.mli) $(MLOBJ_SENS:.cmo=.cmi) 
	$(OCAMLDOC) -g dochtml.cmo \
	    -cvode-doc-root "$(CVODE_DOC_ROOT)" 	\
	    -cvodes-doc-root "$(CVODES_DOC_ROOT)" 	\
	    -ida-doc-root "$(IDA_DOC_ROOT)" 		\
	    -idas-doc-root "$(IDAS_DOC_ROOT)" 		\
	    -kinsol-doc-root "$(KINSOL_DOC_ROOT)" 	\
	    -pp "$(DOCPP)"				\
	    -d ./doc/html/				\
	    -hide Cvode_impl,Ida_impl,Kinsol_impl	\
	    -t "Sundials"				\
	    -intro intro.doc				\
	    $(MLOBJ_MAIN:.cmo=.mli) $(MLOBJ_SENS:.cmo=.mli)

doc/html:
	mkdir $@

# ##

install: sundials.cma sundials.cmxa sundials_wos.cma sundials_wos.cmxa doc META
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

ocamlfind: sundials.cma sundials.cmxa META
	ocamlfind install sundials $(INSTALL_FILES) $(STUBLIBS)

# ##

depend: .depend
.depend:
	$(OCAMLDEP) $(INCLUDES) \
	    -pp "cpp $(CPPFLAGS) -DOCAML_3X=$(OCAML_3X)" \
	    *.mli *.ml > .depend
	$(CC) -MM $(CFLAGS) *.c >> .depend

clean:
	-@(cd examples; make -f Makefile clean)
	-@$(RM) -f $(MLOBJ) $(MLOBJ:.cmo=.cmx) $(MLOBJ:.cmo=.o)
	-@$(RM) -f $(ALL_COBJ) $(MLOBJ:.cmo=.annot)
	-@$(RM) -f $(MLOBJ:.cmo=.cma) $(MLOBJ:.cmo=.cmxa)
	-@$(RM) -f sundials$(XA) sundials_wos$(XA)
	-@$(RM) -f dochtml.cmi dochtml.cmo

cleandoc:
	-@$(RM) -f doc/html/*.html doc/html/style.css

realclean: cleanall
cleanall: clean cleandoc
	-@(cd examples; make -f Makefile cleanall)
	-@$(RM) -f $(MLOBJ:.cmo=.cmi)
	-@$(RM) -f $(MLOBJ:.cmo=.annot)
	-@$(RM) -f sundials.cma sundials.cmxa
	-@$(RM) -f sundials_wos.cma sundials_wos.cmxa
	-@$(RM) -f libmlsundials$(XA) dllmlsundials$(XS)
	-@$(RM) -f libmlsundials_wos$(XA) dllmlsundials_wos$(XS)
	-@$(RM) -f META
	-@$(RM) -f config config.h

-include .depend
