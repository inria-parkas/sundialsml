info:
	@echo "Main targets (cvode, cvodes, ida, idas, kinsol):"
	@echo "  all               compiles all tests without running them"
	@echo "  tests.byte.log    run standard tests (byte code)"
	@echo "  tests.opt.log     run standard tests (native code)"
	@echo "  perf.byte.log     time standard tests (byte code)"
	@echo "  perf.opt.log      time standard tests (native code)"
	@echo "  ocaml             compile the other ocaml examples without running them"

include ../config

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

tests.opt.log tests.byte.log:
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

AWK ?= awk

perf.opt.log perf.byte.log:
	@for s in ${SUBDIRS}; do		\
	    ${MAKE} -C $$s $@;			\
	 done
	@echo ======================================================================
	@grep '^#' $(word 1,$(SUBDIRS))/$@ > $@
	@for s in ${SUBDIRS}; do				\
	    grep '^[^#]' $$s/$@ >> $@;	\
	done >> $@
	cat $@
	@$(if $(findstring $@,$(MAKECMDGOALS)), \
	   echo "Note: \"make $(@:.log=.plot)\" can plot this for you.")

perf.%.dat: perf.%.log
	@sed -e 's#/serial/#/ser/#'		\
	     -e 's#/parallel/#/par/#'	\
	     -e 's#kinsol/#kin/#'	\
	     -e 's#cvodes/#cvs/#'	\
	     -e 's#cvode/#cv/#'	\
	     $< > $@

perf.opt.plot: perf.opt.dat
	TITLE='native code performance' utils/plot.sh $<
	@utils/plot.sh --explain-vars

perf.byte.plot: perf.byte.dat
	TITLE='byte code performance' utils/plot.sh $<
	@utils/plot.sh --explain-vars

perf.opt.%: perf.opt.dat
	TITLE='native code performance' TERMINAL=$(subst perf.opt.,,$@) \
	    OUTPUT=$@ utils/plot.sh $<
	@printf "\nPlot saved in $@.\n"
	@utils/plot.sh --explain-vars

perf.byte.%: perf.byte.dat
	TITLE='byte code performance' TERMINAL=$(subst perf.opt.,,$@) \
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
