include ../../../config

EXAMPLES = kinFerTron_dns.ml   \
	   kinFoodWeb_kry.ml   \
	   kinKrylovDemo_ls.ml \
	   kinLaplace_bnd.ml   \
	   kinRoboKin_dns.ml

ALL_EXAMPLES = $(EXAMPLES)

DIVIDER = "----------------------------------------------------------------------"

##

all: $(EXAMPLES:.ml=.byte)
opt: $(EXAMPLES:.ml=.opt)

lapack: $(LAPACK_EXAMPLES:.ml=.byte)
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

kinFerTron_dns.byte: kinFerTron_dns.ml
kinFerTron_dns.opt:  kinFerTron_dns.ml

kinFoodWeb_kry.byte: kinFoodWeb_kry.ml
kinFoodWeb_kry.opt:  kinFoodWeb_kry.ml

kinKrylovDemo_ls.byte: kinKrylovDemo_ls.ml
kinKrylovDemo_ls.opt:  kinKrylovDemo_ls.ml

kinLaplace_bnd.byte: kinLaplace_bnd.ml
kinLaplace_bnd.opt:  kinLaplace_bnd.ml

kinRoboKin_dns.byte: kinRoboKin_dns.ml
kinRoboKin_dns.opt:  kinRoboKin_dns.ml

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

%.sundials.out : $(EXAMPLESROOT)/kinsol/serial/%
	$< > $@

