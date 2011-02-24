include Makefile.inc

VERSION = 0.5.0

MLOBJ = cvode.cmo 		\
	cvode_nvector.cmo	\
	cvode_serial.cmo 	\
	nvector_array.cmo 	\
	nvector.cmo 		\
	sundials.cmo

COBJ =	cvode_ml$(XO) 		\
	cvode_ml_ba$(XO) 	\
	cvode_ml_nvec$(XO) 	\
	nvector_ml$(XO)

INSTALL_FILES = 		\
    META			\
    $(MLOBJ:.cmo=.cmi)		\
    libmlsundials_cvode$(XA)	\
    sundials_cvode$(XA)		\
    sundials_cvode.cma		\
    sundials_cvode.cmxa

STUBLIBS = dllmlsundials_cvode$(XS)

# ##

.PHONY: all sundials_cvode install doc

all: sundials_cvode.cma sundials_cvode.cmxa

sundials_cvode.cma sundials_cvode.cmxa: $(MLOBJ) $(MLOBJ:.cmo=.cmx) $(COBJ)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS) \
	    -o sundials_cvode -oc mlsundials_cvode $^ \
	    $(LAPACK_LIB) -lsundials_cvode -lsundials_nvecserial

cvode.o: cvode_ml.c
cvode_ml_ba.o: cvode_ml_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(CFLAGS) -DCVODE_ML_BIGARRAYS -o $@ -c $<
cvode_ml_nvec.o: cvode_ml_nvec.c
nvector_ml.o: nvector_ml.c

META: META.in
	@$(ECHO) "version = \"$(VERSION)\"" > $@
	@$(CAT) $< >> $@

doc: doc/html/index.html

doc/html/index.html: $(MLOBJ:.cmo=.mli) $(MLOBJ:.cmo=.cmi) intro.doc
	$(OCAMLDOC) -html		\
	    -pp "$(DOCPP)"		\
	    -d ./doc/html/		\
	    -t "Sundials (CVODE)"	\
	    -intro intro.doc		\
	    $(MLOBJ:.cmo=.mli)

# ##

install: sundials_cvode.cma sundials_cvode.cmxa META
	$(MKDIR) $(PKGDIR)
	$(CP) $(INSTALL_FILES) $(PKGDIR)/
	$(CP) $(STUBLIBS) $(OCAML_INCLUDE)/stublibs/

uninstall:
	for f in $(STUBLIBS); do	 \
	    $(RM) $(OCAML_INCLUDE)/stublibs/$$f || true; \
	done
	for f in $(INSTALL_FILES); do	 \
	    $(RM) $(PKGDIR)/$$f || true; \
	done
	-$(RMDIR) $(PKGDIR)

# ##

depend:
	$(OCAMLDEP) $(INCLUDES) *.mli *.ml > .depend

clean:
	-@$(RM) -f $(MLOBJ) $(MLOBJ:.cmo=.cmx) $(MLOBJ:.cmo=.o)
	-@$(RM) -f $(COBJ) cvode.annot
	-@$(RM) -f $(MLOBJ:.cmo=.cma) $(MLOBJ:.cmo=.cmxa)
	-@$(RM) -f sundials_cvode$(XA)

cleanall: clean
	-@$(RM) -f $(MLOBJ:.cmo=.cmi)
	-@$(RM) -f sundials_cvode.cma sundials_cvode.cmxa
	-@$(RM) -f libmlsundials_cvode$(XA) dllmlsundials_cvode$(XS)
	-@$(RM) -f META

-include .depend
