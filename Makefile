include config

MLOBJ = sundials.cmo nvector.cmo nvector_array.cmo dls.cmo cvode.cmo	\
	spils.cmo spils_nvector.cmo spils_serial.cmo \
	cvode_nvector.cmo cvode_serial.cmo ida.cmo \
	kinsol_nvector.cmo kinsol_serial.cmo \
	ida_nvector.cmo ida_serial.cmo

COMMON_COBJ= sundials_ml$(XO) dls_ml$(XO) nvector_ml$(XO) \

CVODE_COBJ= cvode_ml$(XO) cvode_ml_ba$(XO) cvode_ml_nvec$(XO)

IDA_COBJ= ida_ml$(XO) ida_ml_ba$(XO) ida_ml_nvec$(XO)

KINSOL_COBJ= kinsol_ml$(XO) kinsol_ml_ba$(XO) kinsol_ml_nvec$(XO)

SPILS_COBJ= spils_ml$(XO) spils_ml_ba$(XO) spils_ml_nvec$(XO)

COBJ=$(COMMON_COBJ) $(SPILS_COBJ) $(CVODE_COBJ) $(IDA_COBJ) $(KINSOL_COBJ)

INSTALL_FILES= 			\
    META			\
    $(MLOBJ:.cmo=.cmi)		\
    libmlsundials$(XA)	\
    sundials$(XA)		\
    sundials.cma		\
    sundials.cmxa		\

STUBLIBS=dllmlsundials$(XS)

CFLAGS+=-fPIC

# ##

.PHONY: all sundials install doc

all: sundials.cma sundials.cmxa doc

sundials.cma sundials.cmxa: $(MLOBJ) $(MLOBJ:.cmo=.cmx) $(COBJ)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS) \
	    -o sundials -oc mlsundials $^ \
	    $(OCAML_CVODE_LIBLINK) \
	    $(OCAML_IDA_LIBLINK) \
	    $(OCAML_KINSOL_LIBLINK)

# There are three sets of flags:
#   - one for CVODE-specific files
#   - one for IDA-specific files
#   - one for files common to CVODE and IDA

# The CFLAGS settings for CVODE works for modules common to CVODE and IDA.
$(COMMON_COBJ): %.o: %.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvode_ml.o: cvode_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<
cvode_ml_ba.o: cvode_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) \
	      -DCVODE_ML_BIGARRAYS -o $@ -c $<
cvode_ml_nvec.o: cvode_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

ida_ml.o: ida_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<
ida_ml_ba.o: ida_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) \
	      -DIDA_ML_BIGARRAYS -o $@ -c $<
ida_ml_nvec.o: ida_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

kinsol_ml.o: kinsol_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<
kinsol_ml_ba.o: kinsol_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) \
	      -DKINSOL_ML_BIGARRAYS -o $@ -c $<
kinsol_ml_nvec.o: kinsol_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

spils_ml_ba.o: spils_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) \
	      -DSPILS_ML_BIGARRAYS -o $@ -c $<
spils_ml_nvec.o: spils_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

dochtml.cmo: INCLUDES += -I +ocamldoc
dochtml.cmo: OCAMLFLAGS += -pp "cpp $(CPPFLAGS) -DOCAML_3X=$(OCAML_3X)"

META: META.in
	@$(ECHO) "version = \"$(VERSION)\"" > $@
	@$(CAT) $< >> $@

doc: doc/html/index.html

doc/html/index.html: doc/html dochtml.cmo intro.doc \
		     $(MLOBJ:.cmo=.mli) $(MLOBJ:.cmo=.cmi) 
	$(OCAMLDOC) -g dochtml.cmo \
	    -cvode-doc-root "$(CVODE_DOC_ROOT)" \
	    -ida-doc-root "$(IDA_DOC_ROOT)" \
	    -pp "$(DOCPP)"		\
	    -d ./doc/html/		\
	    -t "Sundials (CVODE, IDA & KINSOL)"	\
	    -intro intro.doc		\
	    $(MLOBJ:.cmo=.mli)

doc/html:
	mkdir $@

# ##

install: sundials.cma sundials.cmxa doc META
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
	$(OCAMLDEP) $(INCLUDES) *.mli *.ml > .depend
	$(CC) -MM $(CFLAGS) *.c >> .depend

clean:
	-@(cd examples; make -f Makefile clean)
	-@$(RM) -f $(MLOBJ) $(MLOBJ:.cmo=.cmx) $(MLOBJ:.cmo=.o)
	-@$(RM) -f $(COBJ) $(MLOBJ:.cmo=.annot)
	-@$(RM) -f $(MLOBJ:.cmo=.cma) $(MLOBJ:.cmo=.cmxa)
	-@$(RM) -f sundials$(XA)
	-@$(RM) -f dochtml.cmi dochtml.cmo

cleandoc:
	-@$(RM) -f doc/html/*.html doc/html/style.css

realclean: cleanall
cleanall: clean
	-@(cd examples; make -f Makefile cleanall)
	-@$(RM) -f $(MLOBJ:.cmo=.cmi)
	-@$(RM) -f sundials.cma sundials.cmxa
	-@$(RM) -f libmlsundials$(XA) dllmlsundials$(XS)
	-@$(RM) -f META

-include .depend
