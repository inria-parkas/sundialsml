
ifeq ($(wildcard config),)
$(error Please run ./configure first)
endif
include config

default: all

# To compile with profiling (with gcc):
# ./configure CFLAGS=-pg OCAMLOPTFLAGS=-p ...

### Common rules
# We can't use pattern rules here because they don't generate
# dependencies and ocamldep doesn't generate dependencies of the form
# foo.cm[iox]: foo.ml or foo.cmi: foo.mli
.SUFFIXES : .mli .ml .cmi .cmo .cmx

.ml.cmo:
	$(OCAMLC) $(OCAMLFLAGS) -c $(INCLUDES) $<

.mli.cmi:
	$(OCAMLC) $(OCAMLFLAGS) -c $(INCLUDES) $<

.ml.cmx:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -c $(INCLUDES) $<

%.o: %.c
	$(CC) -I $(OCAML_INCLUDE) $(CFLAGS) -o $@ -c $<

### Objects shared between sundials.cma and sundials_no_sens.cma.

# Common to CVODE, IDA, KINSOL, and ARKODE.
COBJ_COMMON = sundials_ml$(XO) dls_ml$(XO) $(SLS_ML_XO) nvector_ml$(XO) \
	      spils_ml$(XO) $(NVECPTHREADS_ML_XO) $(NVECOPENMP_ML_XO)

COBJ_MAIN = $(COBJ_COMMON) kinsol_ml$(XO) $(ARKODE_COBJ_MAIN)

MLOBJ_MAIN = sundials_config.cmo sundials.cmo nvector.cmo		\
	     dls_impl.cmo dls.cmo sls_impl.cmo $(SLS_CMO) spils.cmo	\
	     nvector_custom.cmo nvector_array.cmo			\
	     nvector_serial.cmo $(NVECPTHREADS_CMO) $(NVECOPENMP_CMO) 	\
	     cvode_impl.cmo ida_impl.cmo kinsol_impl.cmo		\
	     cvode.cmo kinsol.cmo ida.cmo $(ARKODE_MLOBJ_MAIN)		\
	     $(KLU_MLOBJ_MAIN) $(SUPERLUMT_MLOBJ_MAIN)

CMI_MAIN = $(filter-out sundials_config.cmi,$(filter-out %_impl.cmi,\
	    $(MLOBJ_MAIN:.cmo=.cmi)))

### Objects specific to sundials.cma.
COBJ_SENS  = cvode_ml_s$(XO) ida_ml_s$(XO) cvodes_ml.o idas_ml.o \
	     $(KLU_COBJ_SENS) $(SUPERLUMT_COBJ_SENS)
MLOBJ_SENS = cvodes.cmo idas.cmo \
	     $(KLU_MLOBJ_SENS) $(SUPERLUMT_MLOBJ_SENS)
CMI_SENS = $(MLOBJ_SENS:.cmo=.cmi)

### Objects specific to sundials_no_sens.cma.
COBJ_NO_SENS = cvode_ml$(XO) ida_ml$(XO)			\
	       $(KLU_COBJ_NO_SENS) $(SUPERLUMT_COBJ_NO_SENS)
MLOBJ_NO_SENS =

### Objects specific to sundials_mpi.cma.
COBJ_MPI = nvector_parallel_ml.o kinsol_bbd_ml.o arkode_bbd_ml.o	\
	   cvode_bbd_ml.o cvodes_bbd_ml.o				\
	   ida_bbd_ml.o idas_bbd_ml.o
MLOBJ_MPI = nvector_parallel.cmo kinsol_bbd.cmo arkode_bbd.cmo		\
	    cvode_bbd.cmo cvodes_bbd.cmo				\
	    ida_bbd.cmo idas_bbd.cmo
CMI_MPI = $(MLOBJ_MPI:.cmo=.cmi)

### Other sets of files.

# For `make clean'.  All object files, including ones that may not be
# built/updated under the current configuration.  Duplicates OK.
ALL_COBJ = $(COBJ_MAIN) $(COBJ_SENS) $(COBJ_NO_SENS) $(COBJ_MPI)
ALL_MLOBJ =dochtml.cmo $(MLOBJ_MAIN) $(MLOBJ_SENS) $(MLOBJ_NO_SENS) $(MLOBJ_MPI)
ALL_CMA = sundials.cma sundials_no_sens.cma sundials_mpi.cma sundials_docs.cma

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
	     doc/html/index.html :				\
    INCLUDES += $(MPI_INCLUDES)

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

kinsol_klu_ml.o: kinsol_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

kinsol_superlumt_ml.o: kinsol_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(KINSOL_CFLAGS) -o $@ -c $<

# ARKODE-specific C files.
arkode_ml.o: arkode_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) -o $@ -c $<

arkode_bbd_ml.o: arkode_bbd_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) -o $@ -c $<

arkode_klu_ml.o: arkode_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) -o $@ -c $<

arkode_superlumt_ml.o: arkode_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(ARKODE_CFLAGS) -o $@ -c $<

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

cvode_klu_ml.o: cvode_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvode_klu_ml_s.o: cvode_klu_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) \
	    -o $@ -c $<

cvodes_klu_ml.o: cvodes_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODES_CFLAGS) -o $@ -c $<

cvode_superlumt_ml.o: cvode_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) -o $@ -c $<

cvode_superlumt_ml_s.o: cvode_superlumt_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) $(CVODE_CFLAGS) \
	    -o $@ -c $<

cvodes_superlumt_ml.o: cvodes_superlumt_ml.c
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

ida_klu_ml.o: ida_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

ida_klu_ml_s.o: ida_klu_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) $(IDA_CFLAGS) \
	    -o $@ -c $<

idas_klu_ml.o: idas_klu_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) -o $@ -c $<

ida_superlumt_ml.o: ida_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDA_CFLAGS) -o $@ -c $<

ida_superlumt_ml_s.o: ida_superlumt_ml.c
	$(CC) -DSUNDIALSML_WITHSENS -I $(OCAML_INCLUDE) $(IDA_CFLAGS) \
	    -o $@ -c $<

idas_superlumt_ml.o: idas_superlumt_ml.c
	$(CC) -I $(OCAML_INCLUDE) $(IDAS_CFLAGS) -o $@ -c $<

# Docs.

DOCHTML_PP=$(CPP) $(ML_CPPFLAGS) -DOCAML_3X=$(OCAML_3X)
dochtml.cmo: DOCHTML_PP += -DCVODE_DOC_ROOT=\"$(CVODE_DOC_ROOT_DEFAULT)\"
dochtml.cmo: DOCHTML_PP += -DCVODES_DOC_ROOT=\"$(CVODES_DOC_ROOT_DEFAULT)\"
dochtml.cmo: DOCHTML_PP += -DARKODE_DOC_ROOT=\"$(ARKODE_DOC_ROOT_DEFAULT)\"
dochtml.cmo: DOCHTML_PP += -DIDA_DOC_ROOT=\"$(IDA_DOC_ROOT_DEFAULT)\"
dochtml.cmo: DOCHTML_PP += -DIDAS_DOC_ROOT=\"$(IDAS_DOC_ROOT_DEFAULT)\"
dochtml.cmo: DOCHTML_PP += -DKINSOL_DOC_ROOT=\"$(KINSOL_DOC_ROOT_DEFAULT)\"
dochtml.cmo: DOCHTML_PP += -DMATHJAX_URL=\"$(MATHJAX_URL_DEFAULT)\"
dochtml.cmo: INCLUDES += -I +ocamldoc
dochtml.cmo: OCAMLFLAGS += -pp '$(DOCHTML_PP)'
dochtml.cmo: config
sundials_docs.cma: sundials_config.cmo dochtml.cmo
	$(OCAMLC) $(OCAMLCFLAGS) -o $@ -a $^

META: META.in config
	$(CPP) $(if $(MPI_ENABLED),-DMPI_ENABLED) -DVERSION=\"$(VERSION)\" $< \
	    | grep -v '^#' > $@

doc: doc/html/index.html
	for f in cvode_skel.ml ida_skel.ml kinsol_skel.ml; do \
	    cp examples/ocaml/skeletons/$$f doc/html/; \
	done

DOC_URLS=$(if $(CVODE_DOC_ROOT),-cvode-doc-root "$(CVODE_DOC_ROOT)")	   \
	 $(if $(CVODES_DOC_ROOT),-cvodes-doc-root "$(CVODES_DOC_ROOT)")    \
	 $(if $(ARKODE_DOC_ROOT),-arkode-doc-root "$(ARKODE_DOC_ROOT)")    \
	 $(if $(IDA_DOC_ROOT),-ida-doc-root "$(IDA_DOC_ROOT)")	           \
	 $(if $(IDAS_DOC_ROOT),-idas-doc-root "$(IDAS_DOC_ROOT)")	   \
	 $(if $(KINSOL_DOC_ROOT),-kinsol-doc-root "$(KINSOL_DOC_ROOT)")    \
	 $(if $(MATHJAX_URL),-mathjax "$(MATHJAX_URL)")
doc/html/index.html: OCAML_DOC_ROOT="$(OCAML_DOC_ROOT_DEFAULT)"
doc/html/index.html: doc/html sundials_docs.cma intro.doc		\
		     $(filter-out %_impl.cmi, $(CMI_MAIN))		\
		     $(CMI_SENS) $(if $(MPI_ENABLED), $(CMI_MPI))
	$(OCAMLDOC) -g sundials_docs.cma $(INCLUDES)		\
	    -charset utf-8					\
	    -short-functors					\
	    -colorize-code					\
	    -css-style docstyle.css				\
	    $(DOC_URLS)						\
	    -pp "$(DOCHTML_PP)					\
		-D'OCAML_DOC_ROOT(x)=$(OCAML_DOC_ROOT)x'	\
		-D'VERSION()=$(VERSION)'"			\
	    -d ./doc/html/					\
	    -hide Cvode_impl,Ida_impl,Kinsol_impl,Arkode_impl	\
	    -t "Sundials/ML $(VERSION)p$(VERSIONP)"		\
	    -intro intro.doc					\
	    $(filter-out %_impl.mli, $(CMI_MAIN:.cmi=.mli))	\
	    $(if $(MPI_ENABLED), $(CMI_MPI:.cmi=.mli))		\
	    $(CMI_SENS:.cmi=.mli)

doc/html:
	mkdir $@

doc/html/perf.opt.png: examples/perf.opt.log
	TITLE="OCaml native code performance over C ($(CC) $(CFLAGS))" \
	    TERMINAL=png \
	    FONT=FreeSans,16 \
	    SIZE=2000,1200 \
	    BMARGIN=1600 \
	    DOTSIZE=2.5 DOTTYPE=7 \
	    OUTPUT=$@ examples/utils/plot.sh $<

# Testing the examples

tests.opt.log: examples/tests.opt.log

examples/tests.opt.log: $(INSTALL_CMA:.cma=.cmxa)
	$(MAKE) -C examples PERF_DATA_POINTS=10 $(subst examples/,,$@)

examples/tests.byte.log: $(INSTALL_CMA:.cma=.cmxa)
	$(MAKE) -C examples PERF_DATA_POINTS=10 $(subst examples/,,$@)

examples/perf.opt.log: $(INSTALL_CMA:.cma=.cmxa)
	$(MAKE) -C examples PERF_DATA_POINTS=10 $(subst examples/,,$@)

examples/perf.byte.log: $(INSTALL_CMA:.cma=.cmxa)
	$(MAKE) -C examples PERF_DATA_POINTS=10 $(subst examples/,,$@)

### Install / Uninstall

install: install-sys $(if $(INSTALL_DOCS),install-doc)

# Install to OCaml's system directory -- /usr/lib/ocaml on Debian derivatives.
install-sys: META $(INSTALL_CMA) $(INSTALL_CMA:.cma=.cmxa)
	[ -d $(PKGDIR) ] || $(MKDIR) $(PKGDIR)
	$(CP) $(INSTALL_FILES) $(PKGDIR)
	$(if $(ENABLE_SHARED),[ -d $(STUBDIR) ] || $(MKDIR) $(STUBDIR))
	$(if $(ENABLE_SHARED),$(CP) $(STUBLIBS) $(STUBDIR))

install-doc: doc
	[ -d $(DOCDIR) ] || $(MKDIR) $(DOCDIR)
	[ -d $(DOCDIR)/html ] || $(MKDIR) $(DOCDIR)/html
	$(CP) doc/html/docstyle.css doc/html/*.html $(DOCDIR)/html/

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
	git clone . sandbox
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
	@ # GC bugs that test.opt.log doesn't.
	$(MAKE) -C sandbox/examples PERF_DATA_POINTS=1 perf.opt.log
	$(MAKE) -C sandbox/examples tests.byte.log
	rm -rf sandbox

### Misc

depend: .depend
.depend:
	$(OCAMLDEP) -pp '$(DOCHTML_PP)' *.mli *.ml > .depend
	$(CC) -MM $(CPPFLAGS) *.c >> .depend

clean:
	-@($(MAKE) -C examples clean)
	-@$(RM) -f $(ALL_MLOBJ) $(ALL_MLOBJ:.cmo=.cmx) $(ALL_MLOBJ:.cmo=.o)
	-@$(RM) -f $(ALL_MLOBJ:.cmo=.cmi) $(ALL_MLOBJ:.cmo=.annot) $(ALL_COBJ)
	-@$(RM) -f $(ALL_CMA) $(ALL_CMA:.cma=.cmxa) $(ALL_CMA:.cma=.a)
	-@$(RM) -f $(foreach file,$(INSTALL_CMA:.cma=$(XA)),libml$(file))
	-@$(RM) -f $(foreach file,$(INSTALL_CMA:.cma=$(XS)),dllml$(file))
	-@$(RM) -f $(STUBLIBS)

cleandoc:
	-@$(RM) -f doc/html/*.html

distclean: clean cleandoc
	-@($(MAKE) -C examples distclean)
	-@$(RM) -f META
	-@$(RM) -f config config.h

-include .depend
