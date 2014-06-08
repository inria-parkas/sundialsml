include config

# TODO: compile two .cmas, one with cvodes and one without.

MLOBJ_MAIN = sundials.cmo nvector.cmo nvector_array.cmo dls.cmo \
	     spils.cmo spils_nvector.cmo spils_serial.cmo \
	     cvode.cmo cvode_nvector.cmo cvode_serial.cmo \
	     kinsol.cmo kinsol_nvector.cmo kinsol_serial.cmo \
	     ida.cmo ida_nvector.cmo ida_serial.cmo

MLOBJ_SENS = cvodes_nvector.cmo

MLOBJ_LOCAL = cvode_session_nvector.cmo \
	      cvode_session_serial.cmo

MLOBJ_WOS = $(MLOBJ_MAIN) $(MLOBJ_LOCAL)
MLOBJ = $(MLOBJ_WOS) $(MLOBJ_SENS)

COMMON_COBJ= sundials_ml$(XO) dls_ml$(XO) nvector_ml$(XO) \

CVODE_COBJ= cvode_ml$(XO) cvode_ml_ba$(XO) cvode_ml_nvec$(XO)

CVODES_COBJ= cvodes_ml$(XO) cvodes_ml_ba$(XO) cvodes_ml_nvec$(XO)

IDA_COBJ= ida_ml$(XO) ida_ml_ba$(XO) ida_ml_nvec$(XO)

KINSOL_COBJ= kinsol_ml$(XO) kinsol_ml_ba$(XO) kinsol_ml_nvec$(XO)

SPILS_COBJ= spils_ml$(XO) spils_ml_ba$(XO) spils_ml_nvec$(XO)

COBJ_WOS=$(COMMON_COBJ) $(SPILS_COBJ) $(CVODE_COBJ) $(IDA_COBJ) $(KINSOL_COBJ)
COBJ=$(COBJ_WOS) $(CVODES_COBJ)

INSTALL_FILES= 			\
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

all: sundials.cma sundials.cmxa sundials_wos.cma sundials_wos.cmxa doc

# TODO: fix this:
sundials.cma sundials.cmxa: sundials.cmo sundials.cmx cvode.cmo cvode.cmx \
			    $(MLOBJ_LOCAL) $(MLOBJ_LOCAL:.cmo=.cmx) \
			    $(MLOBJ) $(MLOBJ:.cmo=.cmx) \
			    $(COBJ)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS) \
	    -o sundials -oc mlsundials $^ \
	    $(OCAML_CVODES_LIBLINK) \
	    $(OCAML_IDA_LIBLINK) \
	    $(OCAML_KINSOL_LIBLINK)

# wos = without sensitivity
# TODO: fix this:
sundials_wos.cma sundials_wos.cmxa: sundials.cmo sundials.cmx cvode.cmo cvode.cmx \
				    $(MLOBJ_LOCAL) $(MLOBJ_LOCAL:.cmo=.cmx) \
				    $(MLOBJ_WOS) $(MLOBJ_WOS:.cmo=.cmx) \
				    $(COBJ_WOS)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS) \
	    -o sundials_wos -oc mlsundials_wos $^ \
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

cvode_ml.o: cvode_ml.c spils_ml.h cvode_ml.h sundials_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<
cvode_ml_ba.o: cvode_ml_nvec.c spils_ml.h sundials_ml.h cvode_ml.h nvector_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) \
	      -DCVODE_ML_BIGARRAYS -o $@ -c $<
cvode_ml_nvec.o: cvode_ml_nvec.c spils_ml.h sundials_ml.h cvode_ml.h nvector_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvodes_ml.o: cvodes_ml.c spils_ml.h cvode_ml.h cvodes_ml.h sundials_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) -o $@ -c $<
cvodes_ml_ba.o: cvodes_ml_nvec.c spils_ml.h sundials_ml.h nvector_ml.h \
    		cvode_ml.h cvodes_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) \
	      -DCVODE_ML_BIGARRAYS -o $@ -c $<
cvodes_ml_nvec.o: cvodes_ml_nvec.c spils_ml.h sundials_ml.h nvector_ml.h \
    		  cvode_ml.h cvodes_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

ida_ml.o: ida_ml.c spils_ml.h ida_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<
ida_ml_ba.o: ida_ml_nvec.c nvector_ml.h ida_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) \
	      -DIDA_ML_BIGARRAYS -o $@ -c $<
ida_ml_nvec.o: ida_ml_nvec.c nvector_ml.h ida_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

kinsol_ml.o: kinsol_ml.c spils_ml.h kinsol_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<
kinsol_ml_ba.o: kinsol_ml_nvec.c spils_ml.h nvector_ml.h kinsol_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) \
	      -DKINSOL_ML_BIGARRAYS -o $@ -c $<
kinsol_ml_nvec.o: kinsol_ml_nvec.c spils_ml.h nvector_ml.h kinsol_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

spils_ml_ba.o: spils_ml_nvec.c sundials_ml.h spils_ml.h
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) \
	      -DSPILS_ML_BIGARRAYS -o $@ -c $<
spils_ml_nvec.o: spils_ml_nvec.c sundials_ml.h spils_ml.h
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
	    -cvode-doc-root "$(CVODE_DOC_ROOT)" \
	    -cvodes-doc-root "$(CVODES_DOC_ROOT)" \
	    -ida-doc-root "$(IDA_DOC_ROOT)" \
	    -kinsol-doc-root "$(KINSOL_DOC_ROOT)" \
	    -pp "$(DOCPP)"		\
	    -d ./doc/html/		\
	    -hide Cvode_session_serial,Cvode_session_nvector \
	    -t "Sundials (CVODE, IDA & KINSOL)"	\
	    -intro intro.doc		\
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
	-@$(RM) -f $(COBJ) $(MLOBJ:.cmo=.annot)
	-@$(RM) -f $(MLOBJ:.cmo=.cma) $(MLOBJ:.cmo=.cmxa)
	-@$(RM) -f sundials$(XA) sundials_wos$(XA)
	-@$(RM) -f dochtml.cmi dochtml.cmo

cleandoc:
	-@$(RM) -f doc/html/*.html doc/html/style.css

realclean: cleanall
cleanall: clean
	-@(cd examples; make -f Makefile cleanall)
	-@$(RM) -f $(MLOBJ:.cmo=.cmi)
	-@$(RM) -f sundials.cma sundials.cmxa
	-@$(RM) -f sundials_wos.cma sundials_wos.cmxa
	-@$(RM) -f libmlsundials$(XA) dllmlsundials$(XS)
	-@$(RM) -f libmlsundials_wos$(XA) dllmlsundials_wos$(XS)
	-@$(RM) -f META

-include .depend
