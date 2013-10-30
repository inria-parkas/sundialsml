# Because of the camlp4 trickery, a few things must happen before a .depend can
# be built: the syntax extension (pa_deriving.cmo) must be built; some
# automatically generated .ml files must be created.  Let us call these files
# depend-prereqs.

# generated .ml's; note expr_of_sundials.ml is not generated
GEN_MLS=pprint_sundials.ml pprint_cvode.ml expr_of_cvode.ml \
	pprint_ida.ml expr_of_ida.ml \
	expr_of_cvode_model.ml expr_of_ida_model.ml

# syntax extension
DERIVING_OBJS=camlp4aux.cmo pa_deriving.cmo

DEPEND_PREREQS=$(DERIVING_OBJS) $(GEN_MLS)

# We've tried a few strategies to capture the dependencies of depend-prereqs:
#
# (1) Let the .depend target depend on the depend-prereqs.  When the syntax
#     extension is updated, make tries to update .depend.
#
#     This was unsatisfactory because if the extension is broken (uncompilable)
#     because of an update, then make fails no matter what the command is.
#     This happens even for make clean, which is unacceptable.  We would like
#     make clean to be a reliable command for resetting the whole directory, no
#     matter what state it's in.
#
# (2) Let the .depend target's recipe check if the depend-prereqs are there,
#     and if not, generate a minimal .depend file that covers only the files
#     that can be processed without the extension. Then in the recipes of the
#     syntax extension and generated .ml files, add "$(MAKE) .depend" so that
#     when those files are finally built, the .depend file is automatically
#     updated.
#
#     This didn't work either, because make does not seem to re-read the
#     .depend file when it's updated in the build rules for the syntax
#     extensions or the generated .ml files.  It only seems to re-read .depend
#     when it's updated in the build rule for .depend.  The same make command
#     has to be issued twice to force this to work, but this is
#     counter-intuitive.
#
# The present approach is to split the Makefile into 3 different modes,
# controlled by the variable DEPFILE (or absence thereof).  The Makefile by
# default has only a rule for `clean'.  For any other target, it just
# re-invokes make with the same arguments, but with DEPFILE=.mindepend defined.
# When DEPFILE=.mindepend, this Makefile includes .mindepend, which contains
# dependencies for only the depend-prereqs.  In this mode, the Makefile will
# contain actual build rules for only the depend-prereqs, and for everything
# else re-invokes make with the same arguments but with DEPFILE=.depend.  When
# DEPFILE=.depend, this Makefile enables all build rules and includes .depend
# rather than .mindepend, which lists dependencies for all *.ml and *.mli files
# in this directory.

all: quickcheck_cvode_serial.opt quickcheck_ida_serial.opt
	@echo "### Testing CVODE-serial."
	$(MAKE) -s test-cvode-serial
	@echo ""
	@echo "### Testing IDA-serial."
	$(MAKE) -s test-ida-serial

ifdef DEPFILE

include ../config

SUNDIALS_LIB=-I .. bigarray.cma unix.cma sundials.cma

# Object files that constitute the quickcheck framework, independent of
# sundials.
QUICKCHECK_OBJS=fstream.cmo pprint.cmo quickcheck.cmo

# Object files to be linked into generated CVODE test cases.
CVODE_TEST_OBJS=$(QUICKCHECK_OBJS) pprint_sundials.cmo \
	quickcheck_sundials.cmo pprint_cvode.cmo quickcheck_cvode.cmo

# Object files needed for the generator of CVODE test cases.
CVODE_GEN_OBJS=$(CVODE_TEST_OBJS) camlp4aux.cmo expr_of.cmo \
	expr_of_sundials.cmo expr_of_cvode.cmo expr_of_cvode_model.cmo

# Object files to be linked into generated IDA test cases.
IDA_TEST_OBJS=$(QUICKCHECK_OBJS) pprint_sundials.cmo quickcheck_sundials.cmo \
	pprint_ida.cmo quickcheck_ida.cmo

# Object files needed for the generator of IDA test cases.
IDA_GEN_OBJS=$(IDA_TEST_OBJS) camlp4aux.cmo expr_of.cmo expr_of_sundials.cmo \
	expr_of_ida.cmo expr_of_ida_model.cmo

$(DERIVING_OBJS): %.cmo: %.ml
	$(OCAMLC) $(OCAMLFLAGS) -I +camlp4 -pp $(CAMLP4OF) -o $@ -c $<

$(DERIVING_OBJS:.cmo=.cmx): %.cmx: %.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -I +camlp4 -I .. \
	    -pp $(CAMLP4OF) -o $@ -c $<


ifeq "$(DEPFILE)" ".depend"

test-cvode-serial: quickcheck_cvode_serial.opt $(CVODE_TEST_OBJS)
	./$< $(PARAMS) --exec-file ./test --failed-file failed.ml --compiler \
	    "$(OCAMLC).opt -w -26 -g $(SUNDIALS_LIB) $(filter-out $<,$^)"

quickcheck_cvode_serial.opt: quickcheck_cvode_serial.cmx \
	                     $(CVODE_GEN_OBJS:.cmo=.cmx)

test-ida-serial: quickcheck_ida_serial.opt $(IDA_TEST_OBJS)
	./$< $(PARAMS) --exec-file ./test --failed-file failed.ml --compiler \
	    "$(OCAMLC).opt -w -26 -g $(SUNDIALS_LIB) $(filter-out $<,$^)"

quickcheck_ida_serial.opt: quickcheck_ida_serial.cmx $(IDA_GEN_OBJS:.cmo=.cmx)

P4DERIVING=-pp "$(CAMLP4OF) camlp4aux.cmo pa_deriving.cmo" -I +camlp4

deriving_test.opt: deriving_test.cmx camlp4aux.cmx expr_of.cmx \
	fstream.cmx pprint.cmx

test-deriving: deriving_test.opt
	./$< > deriving_test.out
	diff -u deriving_test.expected deriving_test.out

dump-deriving: pa_deriving.cmo pprint.cmo
	$(CAMLP4OF) camlp4aux.cmo pa_deriving.cmo \
	    -printer Camlp4OCamlPrinter deriving_test.ml

.depend:
	$(OCAMLDEP) $(P4DERIVING) \
	        *.ml *.mli > $@

%.cmo: %.ml
	$(OCAMLC) $(OCAMLFLAGS) $(P4DERIVING) $(SUNDIALS_LIB) -c $< -o $@

%.cmx: %.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) $(P4DERIVING) $(SUNDIALS_LIB:.cma=.cmxa) \
	    -c $< -o $@

%.cmi: %.mli
	$(OCAMLC) $(OCAMLFLAGS) $(P4DERIVING) $(SUNDIALS_LIB) -c $< -o $@

%.opt: %.cmx
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -I +camlp4 \
	    dynlink.cmxa camlp4lib.cmxa \
	    $(SUNDIALS_LIB:.cma=.cmxa) \
	    $(filter-out $<,$^) $< -o $@

-include .depend

else ifeq "$(DEPFILE)" ".mindepend"

.mindepend:
	$(OCAMLDEP) -pp $(CAMLP4OF) $(DERIVING_OBJS:.cmo=.ml) > $@

$(filter *.cmo,$(DEPEND_PREREQS)): %.cmo: %.ml
	$(OCAMLC) -pp $(CAMLP4OF) -c $< -o $@

extract_types.byte: extract_types.cmo
	$(OCAMLC) $(OCAMLFLAGS) -I +camlp4 \
	    dynlink.cma camlp4fulllib.cma camlp4aux.cmo str.cma $^ -o $@

extract_types.cmo: extract_types.ml
	$(OCAMLC) $(OCAMLFLAGS) -pp $(CAMLP4ORF) -I +camlp4 \
	    -I +camlp4/Camlp4Parsers -I +camlp4/Camlp4Printers \
	    -c $< -o $@

SUNDIALS_EXTRACT_EXCLUDE='\.t$$\|array2$$\|val_array$$'

# If you update this target, you should probably update the expr_of_sundials.ml
# target as well.
pprint_sundials.ml: ../sundials.mli extract_types.byte
	./extract_types.byte \
	--header "(* Automatically generated file: do no edit! *)\n\
	(* Pretty-printers for types in the Sundials module. *)\n\
	open Sundials\n\
	open Pprint" \
	--replace "%Tderiving external_types (pretty ~prefix:%M ~opening:%M\
	                ~optional:(Sundials,Roots,RootDirs))" \
	--exclude $(SUNDIALS_EXTRACT_EXCLUDE) \
	$< > $@

IDA_EXTRACT_EXCLUDE="VarTypes.t\\|Id.t\\|root_direction"
# If you update this target, you should probably update the expr_of_ida.ml
# target as well.
pprint_ida.ml: ../ida.mli extract_types.byte
	./extract_types.byte \
	--header "(* Automatically generated file: do no edit! *)\n\
	(* Pretty-printers for types in the Ida module. *)\n\
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

expr_of_ida_model.ml: quickcheck_ida.ml extract_types.byte
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
	Ida.solver_result = solver_result))" \
	--exclude 'IDALang\.' \
	$< > $@


CVODE_EXTRACT_EXCLUDE="root_direction"
# If you update this target, you should probably update the expr_of_cvode.ml
# target as well.
pprint_cvode.ml: ../cvode.mli extract_types.byte
	./extract_types.byte \
	--header "(* Automatically generated file: do no edit! *)\n\
	(* Pretty-printers for types in the Cvode module. *)\n\
	open Cvode\n\
	open Pprint" \
	--replace "%Tderiving external_types (pretty ~prefix:%M ~opening:%M)" \
	--exclude $(CVODE_EXTRACT_EXCLUDE) \
	$< > $@

# If you update this target, you should probably update the pprint_cvode.ml
# target as well.
expr_of_cvode.ml: ../cvode.mli extract_types.byte
	./extract_types.byte \
	--header "(* Automatically generated file: do not edit! *)\n\
	(* expr_of for types defined in the Cvode module. *)\n\
	open Cvode\n\
	open Expr_of\n\
	open Camlp4.PreCast" \
	--replace "%Tderiving external_types (expr_of ~prefix:%M ~opening:%M)" \
	--exclude $(CVODE_EXTRACT_EXCLUDE) \
	$< > $@

expr_of_cvode_model.ml: quickcheck_cvode.ml extract_types.byte
	./extract_types.byte \
	--header "(* Automatically generated file: do not edit! *)\n\
	(* expr_of for types defined in the Quickcheck_cvode module. *)\n\
	open Cvode\n\
	open Expr_of\n\
	open Expr_of_sundials\n\
	open Expr_of_cvode\n\
	open Quickcheck_sundials\n\
	open Quickcheck_cvode\n\
	open Quickcheck\n\
	open Camlp4.PreCast" \
	--replace \
	"%Tderiving external_types \
	(expr_of ~alias:(Carray.t = carray, \
	Cvode.Roots.root_event = root_event, \
	Cvode.root_direction = root_direction, \
	Roots.t = root_info, \
	Cvode.iter = iter, \
	Cvode.lmm = lmm, \
	Cvode.solver_result = solver_result))" \
	--exclude 'CvodeLang\.' \
	$< > $@

-include .mindepend

test-%: $(DEPEND_PREREQS) phony
	$(MAKE) DEPFILE=.depend $@

%.cmo: %.ml $(DEPEND_PREREQS) phony
	$(MAKE) DEPFILE=.depend $@

%.cmi: %.mli $(DEPEND_PREREQS) phony
	$(MAKE) DEPFILE=.depend $@

%.cmx: %.ml $(DEPEND_PREREQS) phony
	$(MAKE) DEPFILE=.depend $@

%.opt: %.ml $(DEPEND_PREREQS) phony
	$(MAKE) DEPFILE=.depend $@

%.byte: %.ml $(DEPEND_PREREQS) phony
	$(MAKE) DEPFILE=.depend $@

else

%:
	@echo "error: Makefile variable DEPFILE is defined but invalid"
	false

endif

else # if DEPFILE is not defined

test-%: phony
	$(MAKE) DEPFILE=.mindepend $@

%.cmo: phony
	$(MAKE) DEPFILE=.mindepend $@

%.cmi: phony
	$(MAKE) DEPFILE=.mindepend $@

%.cmx: phony
	$(MAKE) DEPFILE=.mindepend $@

%.opt: phony
	$(MAKE) DEPFILE=.mindepend $@

%.byte: phony
	$(MAKE) DEPFILE=.mindepend $@

$(GEN_MLS): %.ml: phony
	$(MAKE) DEPFILE=.mindepend $@

endif # ifdef DEPFILE

# generated .ml files should be deleted if the command fails
.DELETE_ON_ERROR:

# clean is always available
clean:
	-rm -f *.cmi *.cmo *.cmx *.o *.opt *.byte *.opt failed test *.out
	-rm -f .depend .mindepend $(GEN_MLS)
	-rm -f phony

phony:
