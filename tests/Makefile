include ../config

test-quickcheck_serial: quickcheck_serial.opt fstream.cmo pprint.cmo \
                        pprint_sundials.cmo quickcheck.cmo
	./$< $(PARAMS) --exec-file ./test --failed-file failed.ml \
	    --compiler "ocamlc -g -I .. bigarray.cma unix.cma \
	                sundials_ida.cma $(filter-out $<,$^)"

test-%: %.opt
	./$< $(PARAMS)

%.cmx: %.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) \
	    -pp camlp4of -I +camlp4 dynlink.cmxa camlp4lib.cmxa \
	    -I .. \
	    bigarray.cmxa unix.cmxa sundials_ida.cmxa -o $@ -c $<
%.cmo: %.ml
	$(OCAMLC) $(OCAMLFLAGS) \
	    -pp camlp4of -I +camlp4 dynlink.cma camlp4lib.cma \
	    -I .. \
	    bigarray.cma unix.cma sundials_ida.cma -o $@ -c $<

%.opt: %.cmx
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    -I +camlp4 dynlink.cmxa camlp4lib.cmxa $(INCLUDES) -I .. \
	    bigarray.cmxa unix.cmxa sundials_ida.cmxa \
	    fstream.cmx pprint.cmx pprint_sundials.cmx quickcheck.cmx $^

.SECONDARY: quickcheck_serial.opt quickcheck.cmx quickcheck.cmo

clean:
	-rm -f *.cmi *.cmo *.cmx *.o *.opt

depend: .depend

.depend:
	$(OCAMLDEP) -pp camlp4of *.ml > $@
-include .depend
