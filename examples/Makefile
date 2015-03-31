info:
	@echo "Main targets (cvode, cvodes, ida, idas, kinsol):"
	@echo "  all               compiles all tests without running them"
	@echo "  tests.byte.log    run standard tests (byte code)"
	@echo "  tests.opt.log     run standard tests (native code)"
	@echo "  perf.byte.log     time standard tests (byte code)"
	@echo "  perf.opt.log      time standard tests (native code)"
	@echo "  ocaml             compile the other ocaml examples without running them"

include ../config

C_EXAMPLES=$(if $(EXAMPLESROOT),,c_examples_unavailable)
c_examples_unavailable:
	@echo "C version of examples not found.  Try running ../configure again"
	@echo "with EXAMPLES=/path/to/sundials/examples.  You can give the examples"
	@echo "directory in the sundials C source tree or (if you installed the C library"
	@echo "with --enable-examples) the one in the installation."
	@false

MODULES  = cvode cvodes ida idas kinsol
VECTYPES = serial $(if $(MPI_ENABLED),parallel)

SUBDIRS=$(foreach m,$(MODULES),$(foreach v,$(VECTYPES),$m/$v))

TOP_TARGETS = tests.byte.log tests.opt.log perf.byte.log perf.opt.log reps
SUBDIR_TARGETS = $(foreach d,$(SUBDIRS),$(foreach t,$(TOP_TARGETS),$d/$t))

.PHONY: default tests.byte.log tests.opt.log tests.byte tests.opt ocaml \
	perf.byte.log perf.opt.log $(SUBDIR_TARGETS)

all: ocaml tests.byte tests.opt

$(SUBDIR_TARGETS):
	@$(MAKE) -C `echo $@ | sed -e 's#/[^/]*$$##'` \
		    `echo $@ | sed -e 's#.*/##'`

tests.opt tests.byte: tests.%: $(foreach s,$(SUBDIR),$s/tests.%)

tests.opt.log tests.byte.log: tests.%.log: $(C_EXAMPLES) \
					$(foreach s,$(SUBDIRS),$s/tests.%.log)
	@echo ======================================================================
	@echo 'Summary (each should be 0):'
	@for s in ${SUBDIRS}; do		\
	    grep "^ *[0-9]* $$s" $$s/$@;	\
	done
	@! for s in ${SUBDIRS}; do		\
	    grep "^ *[0-9]* $$s" $$s/$@;	\
	done | grep -vq '^ *0'

NATIVE_TITLE='OCaml native code performance over C ($(CC) $(CFLAGS))'
BYTE_TITLE  ='OCaml byte code performance over C ($(CC) $(CFLAGS))'

perf.opt.log perf.byte.log: perf.%.log: $(C_EXAMPLES) \
					$(foreach s,$(SUBDIRS),$s/perf.%.log)
	@echo ======================================================================
	@utils/crunchperf -m $(foreach s,$(SUBDIRS),$s/$@) > $@
	utils/crunchperf -s $@
	@$(if $(findstring $@,$(MAKECMDGOALS)), \
	   echo "Note: \"make $(@:.log=.plot)\" can plot this for you.")

reps: $(foreach s,$(SUBDIRS),$s/reps)

NATIVE_TITLE='OCaml native code performance over C ($(CC) $(CFLAGS))'
BYTE_TITLE  ='OCaml byte code performance over C ($(CC) $(CFLAGS))'
PLOTTYPES=jpg png pdf eps
UTILS=utils

perf.opt.plot: perf.opt.log
	TITLE=$(NATIVE_TITLE) utils/plot.sh $<
	@$(UTILS)/plot.sh --explain-vars

perf.byte.plot: perf.byte.log
	TITLE=$(BYTE_TITLE) utils/plot.sh $<
	@$(UTILS)/plot.sh --explain-vars

$(foreach t,$(PLOTTYPES),perf.opt.$t): perf.opt.log
	TITLE=$(NATIVE_TITLE) \
	    TERMINAL=$(subst perf.opt.,,$@)				  \
	    OUTPUT=$@ $(UTILS)/plot.sh $<
	@printf "\nPlot saved in %s.\n" "$@"
	@$(UTILS)/plot.sh --explain-vars

$(foreach t,$(PLOTTYPES),perf.byte.$t): perf.byte.log
	TITLE=$(BYTE_TITLE) \
	    TERMINAL=$(subst perf.byte.,,$@)				  \
	    OUTPUT=$@ $(UTILS)/plot.sh $<
	@printf "\nPlot saved in %s.\n" "$@"
	@$(UTILS)/plot.sh --explain-vars

ocaml:
	${MAKE} -C ocaml

distclean: clean-utils
	@for s in ${SUBDIRS}; do		\
		${MAKE} -C $$s distclean;	\
	done
	@${MAKE} -C ocaml distclean

clean: clean-utils
	@for s in ${SUBDIRS}; do		\
		${MAKE} -C $$s clean;		\
	done
	@${MAKE} -C ocaml clean

clean-utils:
	-@rm -f perf.byte.* perf.opt.*
	-@rm -f $(foreach f,utils/perf utils/crunchperf,\
		    $f $f.cmi $f.cmx $f.cmo $f.o)
