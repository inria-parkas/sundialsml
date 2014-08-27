info:
	@echo "Main targets (cvode, cvodes, ida, idas, kinsol):"
	@echo "  all               compiles all tests without running them"
	@echo "  tests.byte.log    run standard tests (byte code)"
	@echo "  tests.opt.log     run standard tests (native code)"
	@echo "  perf.byte.log     time standard tests, requires GNU utils (byte code)"
	@echo "  perf.opt.log      time standard tests, requires GNU utils (native code)"
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
	@echo "# NB: Examples that took < $(MIN_TIME)s are omitted." > $@
	@echo "# <name>	<OCaml time>	<C time>	<OCaml time / C time>" >> $@
	@for s in ${SUBDIRS}; do						\
	    $(AWK) '/^[^#]/ { if ($$3 >= $(MIN_TIME)) print $$0; }' $$s/$@;	\
	done >> $@
	cat $@
	@$(if $(findstring $@,$(MAKECMDGOALS)), \
	   echo "Note: \"make $(@:.log=.plot)\" can plot this for you.")

perf.%.dat: perf.%.log
	@sed -e 's#/serial/#/serial/\\n#'		\
	     -e 's#/parallel/#/parallel/\\n#'	\
	     -e 's#\.opt##; s#\.byte##'		\
	     $< > $@

BMARGIN=0.1

plot_script =						\
    set bmargin screen $(BMARGIN);			\
    set boxwidth 0.5;					\
    set style fill solid;				\
    set yrange [0:*];					\
    set y2tics 1 nomirror;				\
    set point 2;					\
    plot '$1' using 4:xtic(1) with boxes lc rgb 'red'	\
	title 'OCaml time / C time (left axis)',	\
	'$1' using 3 with points			\
	pointtype 5 lc rgb 'black'			\
	title 'Absolute time for C code (right axis)'	\
	axes x1y2

perf.opt.plot: perf.opt.dat
	@echo "NB: You might have to adjust BMARGIN to fit horizontal axis labels."
	gnuplot -e "set title 'native code performance';\
		    $(call plot_script,$<); pause -1"

perf.byte.plot: perf.byte.dat
	@echo "NB: You might have to adjust BMARGIN to fit horizontal axis labels."
	gnuplot -e "set title 'byte code performance';	\
		    $(call plot_script,$<); pause -1"

perf.opt.%: perf.opt.dat
	gnuplot -e "set title 'native code performance';	\
		    set terminal $(subst perf.opt.,,$@);	\
		    set output '$@';				\
		    $(call plot_script,$<)"
	@echo "Plot saved in $@.  You might have to adjust BMARGIN to"
	@echo "fit horizontal axis labels."

perf.byte.%: perf.byte.dat
	gnuplot -e "set title 'byte code performance';		\
		    set terminal $(subst perf.byte.,,$@);	\
		    set output '$@';				\
		    $(call plot_script,$<)"
	@echo "Plot saved in $@.  You might have to adjust BMARGIN to"
	@echo "fit horizontal axis labels."

ocaml:
	${MAKE} -C ocaml

distclean: clean
clean:
	@for s in ${SUBDIRS}; do		\
		${MAKE} -C $$s clean;		\
	done
	@${MAKE} -C ocaml clean
	-@rm -f perf.byte.* perf.opt.*
