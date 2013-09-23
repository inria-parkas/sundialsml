include config

COMMON_MLOBJ = sundials.cmo        \
	       nvector.cmo         \
	       nvector_array.cmo   \
	       dls.cmo 		
CVODE_MLOBJ= cvode.cmo 		\
	     cvode_nvector.cmo	\
	     cvode_serial.cmo

IDA_MLOBJ= ida.cmo 		\
	    ida_serial.cmo

COMMON_COBJ= sundials_ml.o      \
	     dls_ml.o

CVODE_COBJ= cvode_ml$(XO)       \
	    cvode_ml_ba$(XO) 	\
	    cvode_ml_nvec$(XO) 	\
	    nvector_ml$(XO)

IDA_COBJ= ida_ml$(XO)           \
	  ida_ml_ba$(XO)        \
	  ida_ml_nvec$(XO)      \
	  nvector_ml$(XO)

MLOBJ=$(COMMON_MLOBJ) $(CVODE_MLOBJ) $(IDA_MLOBJ)
COBJ=$(COMMON_COBJ) $(CVODE_COBJ) $(IDA_COBJ)

INSTALL_FILES= 			\
    META			\
    $(MLOBJ:.cmo=.cmi)		\
    libmlsundials_cvode$(XA)	\
    sundials_cvode$(XA)		\
    sundials_cvode.cma		\
    sundials_cvode.cmxa		\
    libmlsundials_ida$(XA)	\
    sundials_ida$(XA)		\
    sundials_ida.cma		\
    sundials_ida.cmxa

STUBLIBS=dllmlsundials_cvode$(XS) dllmlsundials_ida$(XS)

CFLAGS+=-fPIC

# ##

.PHONY: all sundials_cvode install doc

all: sundials_cvode.cma sundials_cvode.cmxa sundials_ida.cma sundials_ida.cmxa

sundials_cvode.cma sundials_cvode.cmxa: $(COMMON_MLOBJ) $(COMMON_MLOBJ:.cmo=.cmx) $(CVODE_MLOBJ) $(CVODE_MLOBJ:.cmo=.cmx) $(COMMON_COBJ) $(CVODE_COBJ)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS) \
	    -o sundials_cvode -oc mlsundials_cvode $^ \
	    $(OCAML_CVODE_LIBLINK)

sundials_ida.cma sundials_ida.cmxa: $(COMMON_MLOBJ) $(COMMON_MLOBJ:.cmo=.cmx) $(IDA_MLOBJ) $(IDA_MLOBJ:.cmo=.cmx) $(COMMON_COBJ) $(IDA_COBJ)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS) \
	    -o sundials_ida -oc mlsundials_ida $^ \
	    $(OCAML_IDA_LIBLINK)

# There three sets of flags:
#   - one for CVODE-specific files
#   - one for IDA-specific files
#   - one for files common to CVODE and IDA
# Is there a way to group these files and specify their flags all at once?

# These modules are common to CVODE and IDA.  The CFLAGS settings for CVODE
# works for these; if it doesn't, there's probably stuff in these files that
# ought to be moved to solver-specific files.
dls_ml.o: dls_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<
sundials_ml.o: sundials_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<
nvector_ml.o: nvector_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvode_ml.o: cvode_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<
cvode_ml_ba.o: cvode_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) \
	      -DCVODE_ML_BIGARRAYS -o $@ -c $<
cvode_ml_nvec.o: cvode_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

ida_nvector.mli:
	$(SED) \
	-e "/^type \(val_array\|der_array\) =/d"			\
	-e "s/ session\( \|\$\)/ 'a session\1/g"			\
	-e "s/\([ (]\)\([^ ]*\) jacobian_arg\([ )]\|\$\)/\1(\2, 'a) jacobian_arg\3/g" \
	-e "s/\([ (]\)val_array\([ )]\|\$\)/\1'a\2/g"			\
	-e "s/\([ (]\)der_array\([ )]\|\$\)/\1'a\2/g"			\
	-e "s/\([ (]\)nvec\([ )]\|\$\)/\1'a nvector\2/g"		\
	-e "s/\([ (]\)solve_arg\([ )]\|\$\)/\1'a solve_arg\2/g"		\
	-e "s/\([ (]\)single_tmp\([ )]\|\$\)/\1'a single_tmp\2/g"	\
	-e "s/\([ (]\)triple_tmp\([ )]\|\$\)/\1'a triple_tmp\2/g"	\
	-e "s/^\(type 'a nvector = \).*/\1'a Nvector.nvector/"		\
	-e "/(\*ENDINTRO\*)/r ida_nvector.doc"			\
	-e "/^(\*STARTINTRO\*)/,/(\*ENDINTRO\*)/d"				\
	$< > $@
ida_ml.o: ida_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<
ida_ml_ba.o: ida_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) \
	      -DIDA_ML_BIGARRAYS -o $@ -c $<
ida_ml_nvec.o: ida_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

dochtml.cmo: INCLUDES += -I +ocamldoc

META: META.in
	@$(ECHO) "version = \"$(VERSION)\"" > $@
	@$(CAT) $< >> $@

doc: doc/html/index.html

doc/html/index.html: doc/html dochtml.cmo \
		     $(MLOBJ:.cmo=.mli) $(MLOBJ:.cmo=.cmi) \
		     intro.doc cvode_nvector.doc
	$(OCAMLDOC) -g dochtml.cmo \
	    -cvode-doc-root "$(CVODE_DOC_ROOT)" \
	    -ida-doc-root "$(IDA_DOC_ROOT)" \
	    -pp "$(DOCPP)"		\
	    -d ./doc/html/		\
	    -t "Sundials (CVODE & IDA)"	\
	    -intro intro.doc		\
	    $(MLOBJ:.cmo=.mli)

doc/html:
	mkdir $@

# ##

install: sundials_cvode.cma sundials_cvode.cmxa META
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

ocamlfind: sundials_cvode.cma sundials_cvode.cmxa META
	ocamlfind install sundials $(INSTALL_FILES) $(STUBLIBS)

# ##

depend: .depend
.depend:
	$(OCAMLDEP) $(INCLUDES) *.mli *.ml > .depend
	$(CC) -MM $(CFLAGS) *.c >> .depend

clean:
	-@(cd examples; make -f Makefile clean)
	-@$(RM) -f $(MLOBJ) $(MLOBJ:.cmo=.cmx) $(MLOBJ:.cmo=.o)
	-@$(RM) -f $(COBJ) cvode.annot
	-@$(RM) -f $(MLOBJ:.cmo=.cma) $(MLOBJ:.cmo=.cmxa)
	-@$(RM) -f sundials_cvode$(XA)
	-@$(RM) -f sundials_ida$(XA)
	-@$(RM) -f dochtml.cmi dochtml.cmo

cleandoc:
	-@$(RM) -f doc/html/*.html doc/html/style.css

realclean: cleanall
cleanall: clean
	-@(cd examples; make -f Makefile cleanall)
	-@$(RM) -f $(MLOBJ:.cmo=.cmi)
	-@$(RM) -f sundials_cvode.cma sundials_cvode.cmxa
	-@$(RM) -f libmlsundials_cvode$(XA) dllmlsundials_cvode$(XS)
	-@$(RM) -f sundials_ida.cma sundials_ida.cmxa
	-@$(RM) -f libmlsundials_ida$(XA) dllmlsundials_ida$(XS)
	-@$(RM) -f META

-include .depend
