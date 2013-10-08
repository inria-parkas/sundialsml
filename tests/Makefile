include ../config

# Object files that need to be preprocessed with camlp4, but which don't use
# the deriving syntax extension.
# extract_types.cmo is excluded from this list because it uses camlp4orf
# rather than camlp4of.
# pa_deriving.cmo is excluded from this list because make .depend must be
# performed after rebuilding pa_deriving.cmo.
CAMLP4_OBJS=camlp4aux.cmo expr_of_sundials.cmo expr_of.cmo

# Object files that implement the deriving extension.  Doesn't include external
# dependencies introduced as a result of preprocessing the deriving extension.
DERIVING_OBJS=camlp4aux.cmo pa_deriving.cmo

# Object files whose sources use the deriving extension.
OBJS_USING_DERIVING=expr_of_ida.cmo expr_of_ida_model.cmo deriving_test.cmo \
	pprint_ida.cmo quickcheck.cmo quickcheck_ida.cmo \
	quickcheck_ida_serial.cmo

CAMLP4_CMIS=expr_of.cmi pprint.cmi

# Object files that constitute the quickcheck framework, independent of
# sundials.
QUICKCHECK_OBJS=fstream.cmo pprint.cmo quickcheck.cmo

# Object files that have to be linked into the test code and the test code
# generator for IDA.
IDA_TEST_OBJS=$(QUICKCHECK_OBJS) quickcheck_sundials.cmo pprint_ida.cmo \
	quickcheck_ida.cmo

# Object files that the IDA test generators need.
IDA_GEN_OBJS=$(IDA_TEST_OBJS) camlp4aux.cmo expr_of.cmo \
	expr_of_sundials.cmo expr_of_ida.cmo expr_of_ida_model.cmo

SUNDIALS_LIB=bigarray.cma unix.cma sundials_ida.cma

test-quickcheck_ida_serial: quickcheck_ida_serial.opt $(IDA_TEST_OBJS)
	./$< $(PARAMS) --exec-file ./test --failed-file failed.ml \
	    --compiler "$(OCAMLC).opt -g -I .. \
	                $(SUNDIALS_LIB) $(filter-out $<,$^)"

quickcheck_ida_serial.byte: quickcheck_ida_serial.cmo $(IDA_GEN_OBJS)

quickcheck_ida_serial.opt: quickcheck_ida_serial.cmx \
                           $(IDA_GEN_OBJS:.cmo=.cmx)

IDA_EXTRACT_EXCLUDE="VarTypes.t\\|Id.t\\|root_direction"
# If you update this target, you should probably update the expr_of_ida.ml
# target as well.
pprint_ida.ml: ../ida.mli extract_types.byte
	./extract_types.byte \
	--header "(* Automatically generated file: do no edit! *)\n\
	(* Pretty-printer for types in the Ida module. *)\n\
	open Ida\n\
	open Pprint" \
	--replace "%Tderiving external_types (pretty ~prefix:%M ~opening:%M)" \
	--exclude $(IDA_EXTRACT_EXCLUDE) \
	$< > $@

# If you update this target, you should probably update the pprint_ida.ml
# target as well.
expr_of_ida.ml: ../ida.mli extract_types.byte
	./extract_types.byte \
	--header "(* Automatically generated file: do not edit! *)\n\
	(* expr_of for types defined in the Ida module. *)\n\
	open Ida\n\
	open Expr_of\n\
	open Camlp4.PreCast" \
	--replace "%Tderiving external_types (expr_of ~prefix:%M ~opening:%M)" \
	--exclude $(IDA_EXTRACT_EXCLUDE) \
	$< > $@

expr_of_ida_model.ml: quickcheck_ida.ml expr_of_sundials.cmo \
	              expr_of_ida.cmo extract_types.byte
	./extract_types.byte \
	--header "(* Automatically generated file: do not edit! *)\n\
	(* expr_of for types defined in the Quickcheck_ida module. *)\n\
	open Ida\n\
	open Expr_of\n\
	open Expr_of_sundials\n\
	open Expr_of_ida\n\
	open Quickcheck_sundials\n\
	open Quickcheck_ida\n\
	open Quickcheck\n\
	open Camlp4.PreCast" \
	--replace \
	"%Tderiving external_types \
	(expr_of ~alias:(Carray.t = carray, \
	Ida.Roots.root_event = root_event, \
	Ida.root_direction = root_direction, \
	Roots.t = root_info, \
	Ida.linear_solver = linear_solver, \
	Ida.solver_result = solver_result))" \
	$< > $@

extract_types.byte: extract_types.cmo
	$(OCAMLC) $(OCAMLFLAGS) -I +camlp4 \
	    dynlink.cma camlp4fulllib.cma camlp4aux.cmo str.cma $^ -o $@

extract_types.cmo: extract_types.ml
	$(OCAMLC) $(OCAMLFLAGS) -pp $(CAMLP4ORF) -I +camlp4 \
	    -I +camlp4/Camlp4Parsers -I +camlp4/Camlp4Printers \
	    -c $< -o $@

test-deriving: deriving_test.opt
	./$< > deriving_test.out
	diff -u deriving_test.expected deriving_test.out

dump-deriving: pa_deriving.cmo pprint.cmo
	$(CAMLP4OF) camlp4aux.cmo pa_deriving.cmo \
	    -printer Camlp4OCamlPrinter \
	    deriving_test.ml

pa_deriving.cmo: pa_deriving.ml
	$(OCAMLC) $(OCAMLFLAGS) -I +camlp4 -I .. -pp $(CAMLP4OF) -o $@ -c $<
	$(MAKE) -B .depend

deriving_test.opt: deriving_test.cmx $(CAMLP4_OBJS:.cmo=.cmx)


test-%: %.opt
	./$< $(PARAMS)

$(CAMLP4_OBJS): %.cmo: %.ml
	$(OCAMLC) $(OCAMLFLAGS) -I +camlp4 -I .. -pp $(CAMLP4OF) -o $@ -c $<

$(CAMLP4_OBJS:.cmo=.cmx): %.cmx: %.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -I +camlp4 -I .. \
	    -pp $(CAMLP4OF) -o $@ -c $<

$(CAMLP4_CMIS): %.cmi: %.mli
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -I +camlp4 -I .. \
	    -pp $(CAMLP4OF) -o $@ -c $<

$(OBJS_USING_DERIVING): %.cmo: %.ml $(DERIVING_OBJS)
	$(OCAMLC) $(OCAMLFLAGS) -I +camlp4 -I .. \
	    -pp "$(CAMLP4OF) $(DERIVING_OBJS)" \
	    $(SUNDIALS_LIB) -o $@ -c $<

$(OBJS_USING_DERIVING:.cmo=.cmx): %.cmx: %.ml $(DERIVING_OBJS)
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -I +camlp4 -I .. \
	    -pp "$(CAMLP4OF) $(DERIVING_OBJS)" \
	    $(SUNDIALS_LIB:.cma=.cmxa) -o $@ -c $<

%.cmo: %.ml
	$(OCAMLC) $(OCAMLFLAGS) -I .. $(SUNDIALS_LIB) -o $@ -c $<

%.cmx: %.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -I .. $(SUNDIALS_LIB:.cma=.cmxa) \
	    -o $@ -c $<

%.opt: %.cmx
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    -I +camlp4 dynlink.cmxa camlp4lib.cmxa $(INCLUDES) -I .. \
	    $(SUNDIALS_LIB:.cma=.cmxa) \
	    $(filter-out $<,$^) $<

%.byte: %.cmo
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    -I +camlp4 dynlink.cma camlp4lib.cma $(INCLUDES) -I .. \
	    $(SUNDIALS_LIB) $(filter-out $<,$^) $<

.SECONDARY: quickcheck_ida_serial.opt quickcheck.cmx quickcheck.cmo \
            $(QUICKCHECK_OBJS)

clean:
	-rm -f *.cmi *.cmo *.cmx *.o *.opt *.byte *.opt failed test *.out

depend: .depend

.depend:
	if test -r camlp4aux.cmo && test -r pa_deriving.cmo; then \
	    $(OCAMLDEP) -pp "$(CAMLP4OF) camlp4aux.cmo pa_deriving.cmo" \
	        *.ml *.mli > $@; \
	else \
	    $(OCAMLDEP) -pp $(CAMLP4OF) $(CAMLP4_OBJS:.cmo=.ml) \
	        pa_deriving.cmo extract_types.cmo \
	        $(CAMLP4_CMIS:.cmi=.mli) > $@; \
	fi

-include .depend

