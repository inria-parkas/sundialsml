include ../config

test-roots: roots.opt
	./$< $(PARAMS)

roots.opt: roots.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../ \
	    bigarray.cmxa unix.cmxa graphics.cmxa \
	    sundials_ida.cmxa $^

clean:
	-rm -f *.cmi *.cmo *.cmx *.o *.opt
