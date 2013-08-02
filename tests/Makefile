include ../config

# Object files pertaining to the quickcheck framework.
QUICKCHECK_OBJS=fstream.cmo pprint.cmo \
                pprint_sundials.cmo quickcheck.cmo quickcheck_ida.cmo

test-quickcheck_serial: quickcheck_serial.opt $(QUICKCHECK_OBJS)
	./$< $(PARAMS) --exec-file ./test --failed-file failed.ml \
	    --compiler "ocamlc -g -I .. bigarray.cma unix.cma \
	                sundials_ida.cma $(filter-out $<,$^)"

quickcheck_serial.opt: $(QUICKCHECK_OBJS:.cmo=.cmx) quickcheck_serial.cmx
failed.opt: $(QUICKCHECK_OBJS:.cmo=.cmx) failed.cmx

test-%: %.opt
	./$< $(PARAMS)

%.cmx: %.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -pp camlp4of -I +camlp4 -I .. \
	    bigarray.cmxa unix.cmxa sundials_ida.cmxa -o $@ -c $<
%.cmo: %.ml
	$(OCAMLC) $(OCAMLFLAGS) -pp camlp4of -I +camlp4 -I .. \
	    bigarray.cma unix.cma sundials_ida.cma -o $@ -c $<

%.opt: %.cmx
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    -I +camlp4 dynlink.cmxa camlp4lib.cmxa $(INCLUDES) -I .. \
	    bigarray.cmxa unix.cmxa sundials_ida.cmxa \
	    $(filter-out $<,$^) $<

.SECONDARY: quickcheck_serial.opt quickcheck.cmx quickcheck.cmo

clean:
	-rm -f *.cmi *.cmo *.cmx *.o *.opt

depend: .depend

.depend:
	$(OCAMLDEP) -pp camlp4of *.ml *.mli > $@

-include .depend
