include ../config

DIRS=cvode cvodes ida idas kinsol
SUBDIRS=serial $(if $(MPI_ENABLED), parallel)

.PHONY: default tests.byte.log tests.opt.log

all:
	@echo "Main targets (cvode, cvodes, ida, idas, kinsol):"
	@echo "  tests.byte.log    run standard tests (byte code)"
	@echo "  tests.opt.log     run standard tests (native code)"

tests.byte.log:
	@for d in ${DIRS}; do		\
	    for s in ${SUBDIRS}; do	\
		${MAKE} -C $$d/$$s $@;	\
	    done;			\
	done
	@echo ======================================================================
	@for d in ${DIRS}; do		\
	    for s in ${SUBDIRS}; do	\
		cat $$d/$$s/$@;		\
	    done;			\
	done

tests.opt.log:
	@for d in ${DIRS}; do		\
	    for s in ${SUBDIRS}; do	\
		${MAKE} -C $$d/$$s $@;	\
	    done;			\
	done
	@echo ======================================================================
	@for d in ${DIRS}; do		\
	    for s in ${SUBDIRS}; do	\
		cat $$d/$$s/$@;		\
	    done;			\
	done

distclean: clean
clean:
	@for d in ${DIRS}; do			\
	    for s in ${SUBDIRS}; do		\
		${MAKE} -C $$d/$$s clean;	\
	    done;				\
	done

