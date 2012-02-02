include ../../config

EXAMPLES = cchatter.byte discontinuous.byte printall.byte \
	   upvseither.byte

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

clean:
	-rm -f $(EXAMPLES:.byte=.cmo) $(EXAMPLES:.byte=.cmx)
	-rm -f $(EXAMPLES:.byte=.o) $(EXAMPLES:.byte=.cmi)

realclean: cleanall
cleanall: clean
	-rm -f $(EXAMPLES) $(EXAMPLES:.byte=.opt)

# #

.SUFFIXES : .ml .byte .opt

.ml.byte:
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    $(INCLUDES) -I ../../ -dllpath ../../ \
	    bigarray.cma unix.cma sundials_cvode.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../ \
	    bigarray.cmxa unix.cmxa sundials_cvode.cmxa $<

