include variables.mk

default: all

# To compile with profiling (with gcc):
# ./configure CFLAGS=-pg OCAMLOPTFLAGS=-p ...

### Build rules.

.PHONY: all sundials install doc clean distclean

all: $(INSTALL_CMA) $(INSTALL_CMA:.cma=.cmxa)

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

sundials_mpi.cma sundials_mpi.cmxa: $(MLOBJ_MPI) $(MLOBJ_MPI:.cmo=.cmx) \
				    $(COBJ_MPI)
	$(OCAMLMKLIB) $(OCAMLMKLIBFLAGS)		\
	    $(if $(ENABLE_SHARED),,-custom)		\
	    -o sundials_mpi -oc mlsundials_mpi $^	\
	    $(LIB_PATH) $(MPI_LIBLINK)

$(MLOBJ_MPI) $(CMI_MPI) $(MLOBJ_MPI:.cmo=.cmx)	\
	     doc/html/index.html :		\
    INCLUDES += $(MPI_INCLUDES)

# The CFLAGS settings for CVODE works for modules common to CVODE and IDA.
$(COBJ_COMMON): %.o: %.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

nvectors/nvector_parallel_ml.o: nvectors/nvector_parallel_ml.c
	$(MPICC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

# KINSOL-specific C files.
kinsol/kinsol_ml.o: kinsol/kinsol_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) $(CSUBDIRS) -o $@ -c $<

kinsol/kinsol_bbd_ml.o: kinsol/kinsol_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) $(CSUBDIRS) -o $@ -c $<

kinsol/kinsol_klu_ml.o: kinsol/kinsol_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) $(CSUBDIRS) -o $@ -c $<

kinsol/kinsol_superlumt_ml.o: kinsol/kinsol_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) $(CSUBDIRS) -o $@ -c $<

# ARKODE-specific C files.
arkode/arkode_ml.o: arkode/arkode_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

arkode/arkode_bbd_ml.o: arkode/arkode_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

arkode/arkode_klu_ml.o: arkode/arkode_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

arkode/arkode_superlumt_ml.o: arkode/arkode_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

# CVODE[S]-specific C files.
cvode/cvode_ml.o: cvode/cvode_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

cvodes/cvode_ml_s.o: cvode/cvode_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(CVODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

cvodes/cvodes_ml.o: cvodes/cvodes_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) \
	    $(CSUBDIRS) -Icvode -o $@ -c $<

cvode/cvode_bbd_ml.o: cvode/cvode_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

cvodes/cvodes_bbd_ml.o: cvodes/cvodes_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) \
	    $(CSUBDIRS) -Icvode -o $@ -c $<

cvode/cvode_klu_ml.o: cvode/cvode_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

cvodes/cvode_klu_ml_s.o: cvode/cvode_klu_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(CVODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

cvodes/cvodes_klu_ml.o: cvodes/cvodes_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) \
	    $(CVODES_CFLAGS) $(CSUBDIRS) -Icvode -o $@ -c $<

cvode/cvode_superlumt_ml.o: cvode/cvode_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

cvodes/cvode_superlumt_ml_s.o: cvode/cvode_superlumt_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(CVODE_CFLAGS) $(CSUBDIRS) -o $@ -c $<

cvodes/cvodes_superlumt_ml.o: cvodes/cvodes_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) \
	    $(CSUBDIRS) -Icvode -o $@ -c $<

# IDA[S]-specific C files.
ida/ida_ml.o: ida/ida_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) $(CSUBDIRS) -o $@ -c $<

idas/ida_ml_s.o: ida/ida_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(IDA_CFLAGS) $(CSUBDIRS) -o $@ -c $<

idas/idas_ml.o: idas/idas_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) $(CSUBDIRS) -Iida -o $@ -c $<

ida/ida_bbd_ml.o: ida/ida_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) $(CSUBDIRS) -o $@ -c $<

idas/idas_bbd_ml.o: idas/idas_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) $(CSUBDIRS) -Iida -o $@ -c $<

ida/ida_klu_ml.o: ida/ida_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) $(CSUBDIRS) -o $@ -c $<

idas/ida_klu_ml_s.o: ida/ida_klu_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(IDA_CFLAGS) $(CSUBDIRS) -o $@ -c $<

idas/idas_klu_ml.o: idas/idas_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) $(CSUBDIRS) -Iida -o $@ -c $<

ida/ida_superlumt_ml.o: ida/ida_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) $(CSUBDIRS) -o $@ -c $<

idas/ida_superlumt_ml_s.o: ida/ida_superlumt_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) \
	    $(IDA_CFLAGS) $(CSUBDIRS) -o $@ -c $<

idas/idas_superlumt_ml.o: idas/idas_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) $(CSUBDIRS) -Iida -o $@ -c $<

# Docs.

### Install / Uninstall

install: install-sys $(if $(INSTALL_DOCS),install-doc)

# Install to OCaml's system directory -- /usr/lib/ocaml on Debian derivatives.
install-sys: META $(INSTALL_CMA) $(INSTALL_CMA:.cma=.cmxa)
	[ -d $(PKGDIR) ] || $(MKDIR) $(PKGDIR)
	$(CP) $(INSTALL_FILES) $(PKGDIR)
	$(if $(ENABLE_SHARED),[ -d $(STUBDIR) ] || $(MKDIR) $(STUBDIR))
	$(if $(ENABLE_SHARED),$(CP) $(STUBLIBS) $(STUBDIR))

install-doc: ../doc
	[ -d $(DOCDIR) ] || $(MKDIR) $(DOCDIR)
	[ -d $(DOCDIR)/html ] || $(MKDIR) $(DOCDIR)/html
	$(CP) ../doc/html/docstyle.css ../doc/html/*.html $(DOCDIR)/html/

install-ocamlfind: install-findlib
install-findlib: META $(INSTALL_CMA) $(INSTALL_CMA:.cma=.cmxa)
	@ocamlfind install sundialsml $(INSTALL_FILES) -optional $(STUBLIBS)


uninstall: uninstall-sys

uninstall-sys:
	-$(RM) -f $(foreach f,$(STUBLIBS),$(STUBDIR)$f)
	-$(RM) -f $(foreach f,$(INSTALL_FILES),$(PKGDIR)$f)
	-$(RMDIR) $(PKGDIR)

uninstall-doc:
	-$(RM) $(DOCDIR)/html/docstyle.css $(DOCDIR)/html/*.html
	-$(RMDIR) $(DOCDIR)/html
	-$(RMDIR) $(DOCDIR)

uninstall-ocamlfind: uninstall-findlib
uninstall-findlib:
	@ocamlfind remove sundialsml

### Preparing Releases

# Sets up a sandbox and makes sure everything compiles.  You need
# everything to do this: gcc, git, OCaml >= 4, OCamlMPI, lapack, etc.
# The sandbox is configured with $(CONFIG_COMMAND), which defaults to
# the command you used to configure the current repository.
distcheck:
	git clone .. sandbox
	cd sandbox && $(CONFIG_COMMAND) > config.log
	@ # Make sure everything is enabled.
	! grep "NOT FOUND" sandbox/config.log
	! grep "DISABLED" sandbox/config.log
	! grep "without lapack" sandbox/config.log
	@ # Tests.  Important ones and short ones first.
	$(MAKE) -C sandbox all
	$(MAKE) -C sandbox/examples tests.opt.log
	$(MAKE) -C sandbox/examples/ocaml
	$(MAKE) -C sandbox doc
	@ # Because perf.opt.log repeats tests, it sometimes uncovers
	@ # GC bugs that tests.opt.log doesn't.
	$(MAKE) -C sandbox/examples PERF_DATA_POINTS=1 perf.opt.log
	$(MAKE) -C sandbox/examples tests.byte.log
	rm -rf sandbox

### Misc

depend: .depend
.depend:
	$(OCAMLDEP) -pp '$(DOCHTML_PP)' *.mli *.ml */*.mli */*.ml > .depend
	$(CC) -MM $(CPPFLAGS) $(SUBDIRS:%=-I %) */*.c >> .depend

clean:
	-@$(RM) -f $(ALL_MLOBJ) $(ALL_MLOBJ:.cmo=.cmx) $(ALL_MLOBJ:.cmo=.o)
	-@$(RM) -f $(ALL_MLOBJ:.cmo=.cmi) $(ALL_MLOBJ:.cmo=.annot) $(ALL_COBJ)
	-@$(RM) -f $(ALL_CMA) $(ALL_CMA:.cma=.cmxa) $(ALL_CMA:.cma=.a)
	-@$(RM) -f $(foreach file,$(INSTALL_CMA:.cma=$(XA)),libml$(file))
	-@$(RM) -f $(foreach file,$(INSTALL_CMA:.cma=$(XS)),dllml$(file))
	-@$(RM) -f $(STUBLIBS)

distclean: clean
	-@$(RM) -f META
	-@$(RM) -f config.h sundials/sundials_config.ml

-include .depend
