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

.PHONY: default tests.byte.log tests.opt.log tests.byte tests.opt ocaml \
	perf.byte.log perf.opt.log

all: ocaml tests.byte tests.opt

tests.opt tests.byte:
	@for s in ${SUBDIRS}; do		\
	    ${MAKE} -C $$s $@;			\
	done

tests.opt.log tests.byte.log: $(C_EXAMPLES)
	@for s in ${SUBDIRS}; do		\
	    ${MAKE} -C $$s $@;			\
	done
	@echo ======================================================================
	@echo 'Summary (each should be 0):'
	@for s in ${SUBDIRS}; do		\
	    grep '^[0-9]' $$s/$@;		\
	done
	@! for s in ${SUBDIRS}; do		\
	    grep '^[0-9]' $$s/$@;		\
	done | grep -q '^[^0]'

# Filter out examples that take less than this to complete.
MIN_TIME ?= 0.2

perf.opt.log perf.byte.log: $(C_EXAMPLES)
	@for s in ${SUBDIRS}; do		\
	    ${MAKE} -C $$s $@;			\
	 done
	@echo ======================================================================
	@for s in ${SUBDIRS}; do				\
	    cat $$s/$@;						\
	    [ $$s = $(lastword $(SUBDIRS)) ] || printf "\n\n";	\
	done > $@
	utils/crunchperf -s $@
	@$(if $(findstring $@,$(MAKECMDGOALS)), \
	   echo "Note: \"make $(@:.log=.plot)\" can plot this for you.")

perf.opt.plot: perf.opt.log
	TITLE='OCaml native code performance over C' utils/plot.sh $<
	@utils/plot.sh --explain-vars

perf.byte.plot: perf.byte.log
	TITLE='OCaml byte code performance over C' utils/plot.sh $<
	@utils/plot.sh --explain-vars

perf.opt.%: perf.opt.log
	TITLE='OCaml native code performance over C' \
	    TERMINAL=$(subst perf.opt.,,$@) \
	    OUTPUT=$@ utils/plot.sh $<
	@printf "\nPlot saved in $@.\n"
	@utils/plot.sh --explain-vars

perf.byte.%: perf.byte.log
	TITLE='OCaml byte code performance over C' \
	    TERMINAL=$(subst perf.opt.,,$@) \
	    OUTPUT=$@ utils/plot.sh $<
	@printf "\nPlot saved in $@.\n"
	@utils/plot.sh --explain-vars

ocaml:
	${MAKE} -C ocaml

distclean: clean
clean:
	@for s in ${SUBDIRS}; do		\
		${MAKE} -C $$s clean;		\
	done
	@${MAKE} -C ocaml clean
	-@rm -f perf.byte.* perf.opt.*
	-@rm -f $(foreach f,utils/perf utils/crunchperf,\
		    $f $f.cmi $f.cmx $f.cmo $f.o)
