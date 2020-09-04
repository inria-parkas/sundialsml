# This file is to be included from the subdirectories.  It is a
# wrapper around examples.mk that adds build rules for
# test_nvector.ml.

EXTRA_DEPS = test_nvector_ml.o test_nvector.cmo

FILES_TO_CLEAN=test_nvector.ml		\
	       test_nvector_ml.c	\
	       test_nvector_ml.o	\
	       test_nvector.o		\
	       test_nvector.annot	\
	       test_nvector.cmi		\
	       test_nvector.cmo		\
	       test_nvector.cmx		\
	       test_nvector.cmt

OCAMLFLAGS += -I .. -custom
OCAMLOPTFLAGS += -I ..

include ../../examples.mk

# use local copies to avoid problems with make -j
test_nvector.ml: ../test_nvector.ml
	cp $< $@; chmod ugo-w $@
test_nvector_ml.c: ../test_nvector_ml.c
	cp $< $@; chmod ugo-w $@

test_nvector.cmo: test_nvector.ml
	$(OCAMLC) $(OCAMLFLAGS) -c $(SUBDIRS:%=-I $(SRC)/%) -o $@ $<

test_nvector.cmx: test_nvector.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -c $(SUBDIRS:%=-I $(SRC)/%) -o $@ $<
