include ../../../config

EXAMPLES = cchatter.byte discontinuous.byte printall.byte \
	   upvseither.byte integr_backward.byte one_over_x.byte \
	   ramp.byte

all: $(EXAMPLES)
opt: $(EXAMPLES:.byte=.opt)

cchatter.byte: cchatter.ml
cchatter.opt: cchatter.ml

discontinuous.byte: discontinuous.ml
discontinuous.opt: discontinuous.ml

printall.byte: printall.ml
printall.opt: printall.ml

sliding.byte: sliding.ml
sliding.opt: sliding.ml

ramp.byte: ramp.ml
ramp.opt: ramp.ml

clean:
	-rm -f $(EXAMPLES:.byte=.cmo) $(EXAMPLES:.byte=.cmx)
	-rm -f $(EXAMPLES:.byte=.o) $(EXAMPLES:.byte=.cmi)

distclean: clean
	-rm -f $(EXAMPLES) $(EXAMPLES:.byte=.opt)

# #

.SUFFIXES : .ml .byte .opt

.ml.byte:
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ -dllpath ../../../ \
	    bigarray.cma unix.cma sundials.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ \
	    bigarray.cmxa unix.cmxa sundials.cmxa $<
