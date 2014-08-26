info:
	@echo "Main targets (cvode, cvodes, ida, idas, kinsol):"
	@echo "  all               compiles all tests without running them"
	@echo "  tests.byte.log    run standard tests (byte code)"
	@echo "  tests.opt.log     run standard tests (native code)"
	@echo "  ocaml             compile the other ocaml examples without running them"

include ../config

MODULES  = cvode cvodes ida idas kinsol
VECTYPES = serial $(if $(MPI_ENABLED),parallel)

SUBDIRS=$(foreach m,$(MODULES),$(foreach v,$(VECTYPES),$m/$v))

.PHONY: default tests.byte.log tests.opt.log tests.byte tests.opt ocaml

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
	    cat $$s/$@				\
	    | grep '^[0-9]';			\
	done
	@! for s in ${SUBDIRS}; do		\
	    cat $$s/$@				\
	    | grep '^[0-9]';			\
	done | grep -q '^[^0]'

ocaml:
	${MAKE} -C ocaml

distclean: clean
clean:
	@for s in ${SUBDIRS}; do		\
		${MAKE} -C $$s clean;		\
	done
	@${MAKE} -C ocaml clean

