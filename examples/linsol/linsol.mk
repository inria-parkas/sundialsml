# This file is to be included from the subdirectories.  It is a
# wrapper around examples.mk that adds build rules for
# test_linsol.ml.

EXTRA_DEPS = test_linsol.cmo

FILES_TO_CLEAN=test_linsol.ml		\
	       test_linsol.o		\
	       test_linsol.annot	\
	       test_linsol.cmi		\
	       test_linsol.cmo		\
	       test_linsol.cmx		\
	       test_linsol.cmt

OCAMLFLAGS += -I .. -custom
OCAMLOPTFLAGS += -I ..

include $(SRCROOT)/examples/examples.mk

# use local copies to avoid problems with make -j
test_linsol.ml: $(SRCROOT)/examples/linsol/test_linsol.ml
	cp $< $@ && chmod ugo-w $@

test_linsol.cmo: test_linsol.ml
	$(OCAMLC) $(OCAMLFLAGS) -c $(SUBDIRS:%=-I $(SRC)/%) -o $@ $<

test_linsol.cmx: test_linsol.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -c $(SUBDIRS:%=-I $(SRC)/%) \
	    $(LIB_PATH:%=-ccopt %) -o $@ $<

