include ../../../config

SRCROOT = ../../../src

all: $(if $(OCAML_GRAPHICS_PATH),dpendulum.byte dpendulum.opt,nographics)

test-pendulum: dpendulum.opt
	./$< $(PARAMS)

nographics:
	@echo "This example requires the OCaml graphics package."
	@echo "(try: opam install graphics, then rerun configure/make.)"

dpendulum.byte: dpendulum.ml
	$(OCAMLC) $(OCAMLFLAGS) -o $@ -dllpath $(SRCROOT) $(INCLUDES) \
	    -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) \
	    -I "$(OCAML_GRAPHICS_PATH)" \
	    unix.cma $(BIGARRAY_CMA) graphics.cma sundials.cma \
	    dpendulum.ml

dpendulum.opt: dpendulum.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) \
	    $(LIB_PATH:%=-ccopt %) -I "$(OCAML_GRAPHICS_PATH)" \
	    unix.cmxa $(BIGARRAY_CMXA) graphics.cmxa sundials.cmxa \
	    dpendulum.ml

clean:
	-@rm -f dpendulum.annot dpendulum.cmt dpendulum.cmti
	-@rm -f dpendulum.cmi dpendulum.cmo dpendulum.cmx dpendulum.o

distclean: clean
	-@rm -f dpendulum.byte dpendulum.opt

