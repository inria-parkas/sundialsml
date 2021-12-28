include ../src/variables.mk

doc: html/index.html
	for f in $(shell ls ../examples/ocaml/skeletons/*_skel.ml); do	\
	    cp $$f html/;	\
	done

DOCHTML_PP=$(CPP) $(ML_CPPFLAGS)
DOCHTML_PLUGIN_PP=$(DOCHTML_PP)						\
                  -DSUNDIALS_DOC_ROOT=\"$(SUNDIALS_DOC_ROOT_DEFAULT)\"	\
                  -DMATHJAX_URL=\"$(MATHJAX_URL_DEFAULT)\"
dochtml.cmo: INCLUDES += -I +ocamldoc -I ../src/sundials
dochtml.cmo: OCAMLFLAGS += -pp '$(DOCHTML_PLUGIN_PP)'
dochtml.cmo: ../config ../src/sundials/sundials_configuration.cmo

dochtml.cmx: INCLUDES += -I +ocamldoc -I ../src/sundials
dochtml.cmx: OCAMLOPTFLAGS += -pp '$(DOCHTML_PLUGIN_PP)'
dochtml.cmx: ../config ../src/sundials/sundials_configuration.cmx

SUNDIALS_DOCS=sundials_docs$(OCAMLDOC_PLUGIN)

sundials_docs.cma: ../src/sundials/sundials_configuration.cmo dochtml.cmo
	$(OCAMLC) $(OCAMLFLAGS) -o $@ -a $^

sundials_docs.cmxs: ../src/sundials/sundials_configuration.cmx dochtml.cmx
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -shared -o $@ $^

../src/sundials/sundials_configuration.cmo ../src/sundials/sundials_configuration.cmx:
	@echo 1>&2 "Error: the binding doesn't seem to be built.  Build the binding"
	@echo 1>&2 "first by typing 'make' at the root of the source tree."
	@false

DOC_URLS=$(if $(SUNDIALS_DOC_ROOT),-sundials-doc-root "$(SUNDIALS_DOC_ROOT)") \
	 $(if $(MATHJAX_URL),-mathjax "$(MATHJAX_URL)")
DOC_SOURCES=$(filter-out %_impl.cmi, $(CMI_MAIN))			\
	    $(CMI_SENS) $(if $(MPI_ENABLED), $(CMI_MPI))		\
	    $(if $(OPENMP_ENABLED),../src/nvectors/nvector_openmp.mli)	\
	    $(if $(PTHREADS_ENABLED),../src/nvectors/nvector_pthreads.mli)
html/index.html: OCAML_DOC_ROOT="$(OCAML_DOC_ROOT_DEFAULT)"
html/index.html: INCLUDES += $(MPI_INCLUDES)
html/index.html: html $(SUNDIALS_DOCS) intro.doc
	$(OCAMLDOC) -html -g $(SUNDIALS_DOCS) $(INCLUDES) $(SUBDIRS:%=-I ../src/%) \
	    -charset utf-8						     \
	    -short-functors						     \
	    -colorize-code						     \
	    -css-style docstyle.css					     \
	    $(DOC_URLS)							     \
	    -pp "$(DOCHTML_PP)						     \
		-D'OCAML_DOC_ROOT(x)=$(OCAML_DOC_ROOT)x'		     \
		-D'VERSION()=$(VERSION)'"				     \
	    -d ./html/							     \
	    -hide Cvode_impl,Ida_impl,Kinsol_impl,Arkode_impl,Sundials_impl,Sundials_LinearSolver_impl,Sundials_NonlinearSolver_impl \
	    -t "Sundials/ML $(VERSION)p$(VERSIONP)"			     \
	    -intro intro.doc						     \
	    $(foreach f, $(DOC_SOURCES:.cmi=.mli), ../src/$f)

html:
	mkdir $@

C_TITLE=C ($(CC) $(filter-out -I% -DNDEBUG=1,$(CFLAGS)))
NATIVE_TITLE=OCaml ($(OCAML_VERSION)) native code performance over $(C_TITLE)
CONFIDENCE=$(shell ../examples/utils/crunchperf --conf ../examples/perf-intv.opt.log)

html/perf.opt.png: ../examples/perf.opt.log
	TITLE="$(NATIVE_TITLE)" \
	    TERMINAL='png noenhanced'		\
	    FONT=FreeSans,16			\
	    SIZE=2000,1200			\
	    BMARGIN=1600			\
	    DOTSIZE=1.5 DOTTYPE=7		\
	    OUTPUT=$@ ../examples/utils/plot.sh $<

html/perf-intv.opt.png: ../examples/perf-intv.opt.log
	TITLE="$(NATIVE_TITLE), $(CONFIDENCE)% confidence" \
	    STYLE=interval			\
	    TERMINAL='png noenhanced'		\
	    FONT=FreeSans,16			\
	    SIZE=2000,1200			\
	    BMARGIN=1600			\
	    DOTSIZE=1.5 DOTTYPE=7		\
	    OUTPUT=$@ ../examples/utils/plot.sh $<

html/perf-intv.opt.pdf: ../examples/perf-intv.opt.log
	TITLE="$(NATIVE_TITLE), $(CONFIDENCE)% confidence" \
	    STYLE=interval			\
	    TERMINAL='pdf noenhanced'		\
	    FONT=FreeSans,16			\
	    SIZE=2000,1200			\
	    BMARGIN=1600			\
	    DOTSIZE=1.5 DOTTYPE=7		\
	    OUTPUT=$@ ../examples/utils/plot.sh $<

install: html/index.html
	[ -d $(DOCDIR) ] || $(MKDIR) $(DOCDIR)
	[ -d $(DOCDIR)html ] || $(MKDIR) $(DOCDIR)html
	$(CP) ../doc/html/docstyle.css ../doc/html/*.html $(DOCDIR)html/

uninstall:
	-$(RM) $(DOCDIR)html/docstyle.css $(DOCDIR)html/*.html
	-$(RMDIR) $(DOCDIR)html
	-$(RMDIR) $(DOCDIR)

distclean clean:
	-@$(RM) -f html/*.html sundials_docs.cmxs sundials_docs.cma
	-@$(RM) -f dochtml.cmo dochtml.cmx dochtml.cmi dochtml.cmt dochtml.o

depend: .depend
.depend:
	$(OCAMLDEP) -pp '$(DOCHTML_PP)' $(SUBDIRS:%=-I ../src/%) > .depend

-include .depend
