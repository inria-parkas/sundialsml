include ../../config

DIRS=skeletons ball linear misc pendulum/ida sincos

.PHONY: default tests.byte.log tests.opt.log

all:
	@for d in ${DIRS}; do		\
	    ${MAKE} -C $$d $@;		\
	done

cleanall distclean:
	@for d in ${DIRS}; do		\
	    ${MAKE} -C $$d $@;		\
	done

clean:
	@for d in ${DIRS}; do		\
	    ${MAKE} -C $$d $@;		\
	done

