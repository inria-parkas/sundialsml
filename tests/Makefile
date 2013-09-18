include ../config

# Object files that implement the deriving extension.
CAMLP4_OBJS=fstream.cmo camlp4aux.cmo pa_deriving.cmo expr_of.cmo pprint.cmo
CAMLP4_CMIS=expr_of.cmi pprint.cmi

# Object files pertaining to the quickcheck framework.
QUICKCHECK_OBJS=fstream.cmo pprint.cmo camlp4aux.cmo expr_of.cmo \
                quickcheck.cmo quickcheck_sundials.cmo quickcheck_ida.cmo

test-quickcheck_ida_serial: quickcheck_ida_serial.opt $(QUICKCHECK_OBJS)
	./$< $(PARAMS) --exec-file ./test --failed-file failed.ml \
	    --compiler "$(OCAMLC).opt -g -I .. -I +camlp4 \
	        dynlink.cma camlp4lib.cma \
	        bigarray.cma unix.cma sundials_ida.cma $(filter-out $<,$^)"

quickcheck_ida_serial.byte: quickcheck_ida_serial.cmo $(QUICKCHECK_OBJS)

quickcheck_ida_serial.opt: quickcheck_ida_serial.cmx \
                           $(QUICKCHECK_OBJS:.cmo=.cmx)

quickcheck_ida.opt: $(QUICKCHECK_OBJS:.cmo=.cmx) quickcheck_ida.cmx

failed.opt: $(QUICKCHECK_OBJS:.cmo=.cmx) failed.cmx

sundials_types.ml: ../sundials.mli extract_types.byte
	./extract_types.byte $< > $@

extract_types.byte: extract_types.cmo
	$(OCAMLC) $(OCAMLFLAGS) -I +camlp4 \
	    dynlink.cma camlp4fulllib.cma camlp4aux.cmo str.cma $^ -o $@

extract_types.cmo: extract_types.ml
	$(OCAMLC) $(OCAMLFLAGS) -pp $(CAMLP4ORF) -I +camlp4 -I +camlp4/Camlp4Parsers -I +camlp4/Camlp4Printers \
	    -c $< -o $@

deriving.cmo: deriving.ml camlp4aux.cmo
	$(OCAMLC) $(OCAMLFLAGS) -pp $(CAMLP4ORF) -I +camlp4 -c $< -o $@

# The deriving preprocessor depends on these files, so they must be built
# without -pp "camlp4 pa_deriving.cmo"
$(CAMLP4_OBJS): %.cmo: %.ml
	$(OCAMLC) $(OCAMLFLAGS) -c -I +camlp4 -pp $(CAMLP4OF) $< -o $@
$(CAMLP4_CMIS): %.cmi: %.mli
	$(OCAMLC) $(OCAMLFLAGS) -c -I +camlp4 -pp $(CAMLP4OF) $< -o $@

test-deriving: deriving_test.opt
	./$< > deriving_test.out
	diff -u deriving_test.expected deriving_test.out

dump-deriving: pa_deriving.cmo pprint.cmo
	$(CAMLP4OF) camlp4aux.cmo pa_deriving.cmo \
	    -printer Camlp4OCamlPrinter \
	    deriving_test.ml

deriving_test.opt: deriving_test.cmx $(CAMLP4_OBJS:.cmo=.cmx)


test-%: %.opt
	./$< $(PARAMS)

%.cmx: %.ml $(CAMLP4_OBJS)
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -I +camlp4 -I .. \
	    -pp "$(CAMLP4OF) camlp4aux.cmo pa_deriving.cmo" \
	    bigarray.cmxa unix.cmxa sundials_ida.cmxa -o $@ -c $<
%.cmo: %.ml $(CAMLP4_OBJS)
	$(OCAMLC) $(OCAMLFLAGS) -I +camlp4 -I .. \
	    -pp "$(CAMLP4OF) camlp4aux.cmo pa_deriving.cmo" \
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
	-rm -f *.cmi *.cmo *.cmx *.o *.opt *.byte *.opt failed test *.out

depend: .depend

.depend:
	$(OCAMLDEP) -pp $(CAMLP4OF) $(CAMLP4_OBJS:.cmo=.ml) \
	    $(CAMLP4_CMIS:.cmi=.mli) > $@
	$(OCAMLDEP) -pp "$(CAMLP4OF) camlp4aux.cmo pa_deriving.cmo" \
	    $(filter-out $(CAMLP4_OBJS:.cmo=.ml),*.ml) \
	    $(filter-out $(CAMLP4_OBJS:.cmo=.mli),*.mli) >> $@

include .depend

