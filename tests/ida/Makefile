include ../../config

test-%: %.opt
	./$< $(PARAMS)

%.opt: %.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../ \
	    bigarray.cmxa unix.cmxa graphics.cmxa \
	    sundials_ida.cmxa $^

.SECONDARY: experiment.opt roots.opt

clean:
	-rm -f *.cmi *.cmo *.cmx *.o *.opt
