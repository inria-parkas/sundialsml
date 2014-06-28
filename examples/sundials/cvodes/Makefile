include ../../../config

EXAMPLES = cvsAdvDiff_ASAi_bnd.ml \
	   cvsAdvDiff_FSA_non.ml  \
	   cvsDiurnal_FSA_kry.ml  \
	   cvsFoodWeb_ASAi_kry.ml \
	   cvsFoodWeb_ASAp_kry.ml \
	   cvsHessian_ASA_FSA.ml  \
	   cvsRoberts_ASAi_dns.ml \
	   cvsRoberts_FSA_dns.ml

LAPACK_EXAMPLES =

ALL_EXAMPLES = $(EXAMPLES) $(LAPACK_EXAMPLES)

DIVIDER = "----------------------------------------------------------------------"

##

all: byte opt
byte: $(EXAMPLES:.ml=.byte)
opt: $(EXAMPLES:.ml=.opt)

lapack: lapack.byte lapack.opt
lapack.byte: $(LAPACK_EXAMPLES:.ml=.byte)
lapack.opt: $(LAPACK_EXAMPLES:.ml=.opt)

genlog = \
    @(for f in $(1); do  \
	echo $(DIVIDER);                  \
	echo "--$$f";                    \
	cat $$f;                         \
     done;                               \
     echo "Summary (each should be 0):"; \
     for f in $(1); do  \
	wc -l $$f;                       \
     done;                               \
    ) > $(2)

# requires a valid EXAMPLESROOT to avoid:
#   make: *** No rule to make target `cvAdvDiff_bnd.sundials.out',
#					needed by `tests.byte.log'.  Stop.
# e.g., make EXAMPLESROOT=/path/to/sundials-2.x.0/examples tests.log
tests.log: tests.byte.log tests.opt.log

tests.byte.log: $(EXAMPLES:.ml=.byte.out) $(EXAMPLES:.ml=.sundials.out)
	@for f in $(EXAMPLES:.ml=); do \
	    diff -u $$f.sundials.out $$f.byte.out > $$f.byte.diff || true; \
	done
	$(call genlog, $(EXAMPLES:.ml=.byte.diff), $@)
	cat $@

tests.opt.log: $(EXAMPLES:.ml=.opt.out) $(EXAMPLES:.ml=.sundials.out)
	@for f in $(EXAMPLES:.ml=); do \
	    diff -u $$f.sundials.out $$f.opt.out > $$f.opt.diff || true; \
	done
	$(call genlog, $(EXAMPLES:.ml=.opt.diff), $@)
	cat $@

tests.self.log: $(EXAMPLES:.ml=.byte.out) $(EXAMPLES:.ml=.opt.out)
	@for f in $(EXAMPLES:.ml=); do \
	    diff -u $$f.byte.out $$f.opt.out > $$f.self.diff || true; \
	done
	$(call genlog, $(EXAMPLES:.ml=.self.diff), $@)
	cat $@

lapack-tests.byte.log: $(LAPACK_EXAMPLES:.ml=.byte.out) \
		       $(LAPACK_EXAMPLES:.ml=.sundials.out)
	@for f in $(LAPACK_EXAMPLES:.ml=); do \
	    diff -u $$f.sundials.out $$f.byte.out > $$f.byte.diff || true; \
	done
	$(call genlog, $(LAPACK_EXAMPLES:.ml=.byte.diff), $@)
	cat $@

lapack-tests.opt.log: $(LAPACK_EXAMPLES:.ml=.opt.out) \
		      $(LAPACK_EXAMPLES:.ml=.sundials.out)
	@for f in $(LAPACK_EXAMPLES:.ml=); do \
	    diff -u $$f.sundials.out $$f.opt.out > $$f.opt.diff || true; \
	done
	$(call genlog, $(LAPACK_EXAMPLES:.ml=.opt.diff), $@)
	cat $@

##

# Tests with command-line arguments (make same choice as Sundials distrib)

cvsAdvDiff_FSA_non.byte.out: cvsAdvDiff_FSA_non.byte
	./$< -sensi sim t > $@

cvsDiurnal_FSA_kry.byte.out: cvsDiurnal_FSA_kry.byte
	./$< -sensi sim t > $@

cvsRoberts_FSA_dns.byte.out: cvsRoberts_FSA_dns.byte
	./$< -sensi sim t > $@

cvsAdvDiff_FSA_non.opt.out: cvsAdvDiff_FSA_non.opt
	./$< -sensi sim t > $@

cvsDiurnal_FSA_kry.opt.out: cvsDiurnal_FSA_kry.opt
	./$< -sensi sim t > $@

cvsRoberts_FSA_dns.opt.out: cvsRoberts_FSA_dns.opt
	./$< -sensi sim t > $@

##

clean:
	-@rm -f $(ALL_EXAMPLES:.ml=.cmo) $(ALL_EXAMPLES:.ml=.cmx)
	-@rm -f $(ALL_EXAMPLES:.ml=.o) $(ALL_EXAMPLES:.ml=.cmi)
	-@rm -f $(ALL_EXAMPLES:.ml=.c.log) $(ALL_EXAMPLES:.ml=.ml.log)

realclean: cleanall
cleanall: clean
	-@rm -f $(ALL_EXAMPLES:.ml=.byte) $(ALL_EXAMPLES:.ml=.opt)
	-@rm -f $(ALL_EXAMPLES:.ml=.byte.out) $(ALL_EXAMPLES:.ml=.opt.out)
	-@rm -f $(ALL_EXAMPLES:.ml=.sundials.out)
	-@rm -f $(ALL_EXAMPLES:.ml=.annot)
	-@rm -f $(ALL_EXAMPLES:.ml=.byte.diff) $(ALL_EXAMPLES:.ml=.opt.diff)
	-@rm -f $(ALL_EXAMPLES:.ml=.self.diff)
	-@rm -f tests.log lapack-tests.log tests.self.log
	-@rm -f tests.byte.log lapack-tests.byte.log
	-@rm -f tests.opt.log lapack-tests.opt.log

# #

.SUFFIXES : .ml .byte .opt .diff .sundials.out

.ml.byte:
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ -dllpath ../../../ bigarray.cma unix.cma \
	    sundials.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ bigarray.cmxa unix.cmxa \
	    sundials.cmxa $<

%.byte.out: %.byte
	./$< > $@

%.opt.out: %.opt
	./$< > $@

%.sundials.out : $(EXAMPLESROOT)/cvodes/serial/%
	$< > $@

