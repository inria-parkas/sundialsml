
DIRS=`find . -mindepth 1 -maxdepth 1 -type d \! -name .svn`

.PHONY: all clean realclean

all:
	@for i in ${DIRS}; do ${MAKE} -C $$i all; done;
clean:
	@for i in ${DIRS}; do ${MAKE} -C $$i clean; done;
realclean cleanall: clean
	@for i in ${DIRS}; do ${MAKE} -C $$i realclean; done;
