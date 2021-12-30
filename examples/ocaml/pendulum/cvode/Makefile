include ../../../../config

SRCROOT = ../../../../src

all: $(if $(OCAML_GRAPHICS_PATH),pendulum.byte pendulum.opt,nographics)

test-pendulum: pendulum.opt
	./$< $(PARAMS)

nographics:
	@echo "This example requires the OCaml graphics package."
	@echo "(try: opam install graphics, then rerun configure/make.)"

pendulum.byte: showpendulum.ml pendulum.ml
	$(OCAMLC) $(OCAMLFLAGS) -o $@ -dllpath $(SRCROOT) $(INCLUDES) \
	    -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) \
	    -I "$(OCAML_GRAPHICS_PATH)" \
	    unix.cma bigarray.cma graphics.cma sundials.cma \
	    showpendulum.ml pendulum.ml

pendulum.opt: showpendulum.cmx pendulum.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) \
	    $(LIB_PATH:%=-ccopt %) -I "$(OCAML_GRAPHICS_PATH)" \
	    unix.cmxa bigarray.cmxa graphics.cmxa sundials.cmxa \
	    showpendulum.ml pendulum.ml

showpendulum.cmo: showpendulum.ml
showpendulum.cmx: showpendulum.ml

clean:
	-@rm -f showpendulum.cmi showpendulum.cmo showpendulum.cmx showpendulum.o
	-@rm -f showpendulum.annot pendulum.annot
	-@rm -f showpendulum.cmt showpendulum.cmti
	-@rm -f pendulum.cmi pendulum.cmo pendulum.cmx pendulum.o
	-@rm -f pendulum.cmt pendulum.cmti

distclean: clean
	-@rm -f pendulum.byte pendulum.opt

