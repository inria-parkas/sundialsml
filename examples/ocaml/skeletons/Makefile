include ../../../config

# TODO: finish the other skeletons and add them to the docs.

SRCROOT = ../../..

EXAMPLES = cvode_skel.byte  \
	   ida_skel.byte    \
	   arkode_skel.byte \
	   kinsol_skel.byte

all: $(EXAMPLES)
opt: $(EXAMPLES:.byte=.opt)

cvode_skel.byte: cvode_skel.ml
cvode_skel.opt: cvode_skel.ml

ida_skel.byte: ida_skel.ml
ida_skel.opt: ida_skel.ml

kinsol_skel.byte: kinsol_skel.ml
kinsol_skel.opt: kinsol_skel.ml

intro1.out: intro1.ml
	# Requires that ../../.. is called "sundialsml"
	OCAMLPATH=../../../.. $(OCAML) -stdin < $^ > $@

intro2.out: intro2.ml
	# Requires that sundialsml is installed
	ocamlfind ocamlmktop -o partop -package sundialsml.mpi,findlib -linkpkg
	OCAMLPATH=../../../.. $(OCAML) -I $(OCAMLMPI) mpi.cma -stdin < $^ > $@

# #

clean:
	-@rm -f $(EXAMPLES:.byte=.cmo) $(EXAMPLES:.byte=.cmx)
	-@rm -f $(EXAMPLES:.byte=.o) $(EXAMPLES:.byte=.cmi)
	-@rm -f $(EXAMPLES:.byte=.annot) $(EXAMPLES)
	-@rm -f intro1.out intro2.out

distclean: clean
	-@rm -f $(EXAMPLES) $(EXAMPLES:.byte=.opt)

# #

.SUFFIXES : .ml .byte .opt

.ml.byte:
	$(OCAMLC) $(OCAMLFLAGS) -o $@ -dllpath $(SRCROOT) \
	    $(INCLUDES) -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) \
	    bigarray.cma sundials.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ -dllpath $(SRCROOT) \
	    $(INCLUDES) -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) \
	    bigarray.cmxa sundials.cmxa $<

