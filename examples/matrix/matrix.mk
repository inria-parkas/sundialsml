# This file is to be included from the subdirectories.  It is a
# wrapper around examples.mk that adds build rules for
# test_matrix.ml.

EXTRA_DEPS = test_matrix.cmo

FILES_TO_CLEAN=test_matrix.ml		\
	       test_matrix.o		\
	       test_matrix.annot	\
	       test_matrix.cmi		\
	       test_matrix.cmo		\
	       test_matrix.cmx		\
	       test_matrix.cmt

OCAMLFLAGS += -I .. -custom
OCAMLOPTFLAGS += -I ..

include ../../examples.mk

# use local copies to avoid problems with make -j
test_matrix.ml: ../test_matrix.ml
	cp $< $@ && chmod ugo-w $@

test_matrix.cmo: test_matrix.ml
	$(OCAMLC) $(OCAMLFLAGS) -c $(SUBDIRS:%=-I $(SRC)/%) -o $@ $<

test_matrix.cmx: test_matrix.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -c $(SUBDIRS:%=-I $(SRC)/%) \
	    $(LIB_PATH:%=-ccopt %) -o $@ $<

