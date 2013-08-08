include ../config

# Object files pertaining to the quickcheck framework.
QUICKCHECK_OBJS=fstream.cmo pprint.cmo \
                quickcheck.cmo quickcheck_sundials.cmo quickcheck_ida.cmo

test-quickcheck_ida_serial: quickcheck_ida_serial.opt $(QUICKCHECK_OBJS)
	./$< $(PARAMS) --exec-file ./test --failed-file failed.ml \
	    --compiler "ocamlc -g -I .. bigarray.cma unix.cma \
	                sundials_ida.cma $(filter-out $<,$^)"

quickcheck_ida_serial.byte: quickcheck_ida_serial.cmo $(QUICKCHECK_OBJS)

quickcheck_ida_serial.opt: quickcheck_ida_serial.cmx \
                           $(QUICKCHECK_OBJS:.cmo=.cmx)

quickcheck_ida.opt: $(QUICKCHECK_OBJS:.cmo=.cmx) quickcheck_ida.cmx

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

%.byte: %.cmo
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    -I +camlp4 dynlink.cma camlp4lib.cma $(INCLUDES) -I .. \
	    bigarray.cma unix.cma sundials_ida.cma \
	    $(filter-out $<,$^) $<

.SECONDARY: quickcheck_ida_serial.opt quickcheck.cmx quickcheck.cmo \
            $(QUICKCHECK_OBJS)

clean:
	-rm -f *.cmi *.cmo *.cmx *.o *.opt

depend: .depend

.depend:
	$(OCAMLDEP) -pp camlp4of *.ml *.mli > $@

-include .depend
