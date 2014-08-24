include ../../../config

EXAMPLES = kinFerTron_dns.ml   		\
	   kinFoodWeb_kry.ml   		\
	   kinFoodWeb_kry_custom.ml     \
	   kinKrylovDemo_ls.ml 		\
	   kinLaplace_bnd.ml   		\
	   kinRoboKin_dns.ml		\
	   kinRoboKin_dns_alt.ml 	\

ALL_EXAMPLES = $(EXAMPLES) kinFoodWeb_kry_custom.ml

DIVIDER = "----------------------------------------------------------------------"

##

all: byte opt

byte: $(ALL_EXAMPLES:.ml=.byte)
opt: $(ALL_EXAMPLES:.ml=.opt)

genlog = \
    @(for f in $(1); do  			\
	echo $(DIVIDER);                 	\
	echo "--kinsol/serial/$$f";      	\
	cat $$f;                         	\
     done;                               	\
     echo "Summary (each should be 0):"; 	\
     for f in $(1); do  			\
	(cd ../..; wc -l kinsol/serial/$$f);	\
     done;                               	\
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
	@! grep '^[0-9]' $@ | grep -q '^[^0]'

tests.opt.log: $(EXAMPLES:.ml=.opt.out) $(EXAMPLES:.ml=.sundials.out)
	@for f in $(EXAMPLES:.ml=); do \
	    diff -u $$f.sundials.out $$f.opt.out > $$f.opt.diff || true; \
	done
	$(call genlog, $(EXAMPLES:.ml=.opt.diff), $@)
	cat $@
	@! grep '^[0-9]' $@ | grep -q '^[^0]'

tests.self.log: $(EXAMPLES:.ml=.byte.out) $(EXAMPLES:.ml=.opt.out)
	@for f in $(EXAMPLES:.ml=); do \
	    diff -u $$f.byte.out $$f.opt.out > $$f.self.diff || true; \
	done
	$(call genlog, $(EXAMPLES:.ml=.self.diff), $@)
	cat $@
	@! grep '^[0-9]' $@ | grep -q '^[^0]'

##

kinFerTron_dns.byte: kinFerTron_dns.ml
kinFerTron_dns.opt:  kinFerTron_dns.ml

kinFoodWeb_kry.byte: kinFoodWeb_kry.ml
kinFoodWeb_kry.opt:  kinFoodWeb_kry.ml

kinFoodWeb_kry_custom.byte: kinFoodWeb_kry_custom.ml
kinFoodWeb_kry_custom.opt:  kinFoodWeb_kry_custom.ml
kinFoodWeb_kry_custom.sundials.out: kinFoodWeb_kry.sundials.out
	cp $< $@

kinKrylovDemo_ls.byte: kinKrylovDemo_ls.ml
kinKrylovDemo_ls.opt:  kinKrylovDemo_ls.ml

kinLaplace_bnd.byte: kinLaplace_bnd.ml
kinLaplace_bnd.opt:  kinLaplace_bnd.ml

kinRoboKin_dns.byte: kinRoboKin_dns.ml
kinRoboKin_dns.opt:  kinRoboKin_dns.ml

kinRoboKin_dns_alt.byte: kinRoboKin_dns_alt.ml
kinRoboKin_dns_alt.opt:  kinRoboKin_dns_alt.ml
kinRoboKin_dns_alt.sundials.out: kinRoboKin_dns.sundials.out
	cp $< $@

##

distclean: clean
clean:
	-@rm -f $(ALL_EXAMPLES:.ml=.cmo) $(ALL_EXAMPLES:.ml=.cmx)
	-@rm -f $(ALL_EXAMPLES:.ml=.o) $(ALL_EXAMPLES:.ml=.cmi)
	-@rm -f $(ALL_EXAMPLES:.ml=.c.log) $(ALL_EXAMPLES:.ml=.ml.log)
	-@rm -f $(ALL_EXAMPLES:.ml=.byte) $(ALL_EXAMPLES:.ml=.opt)
	-@rm -f $(ALL_EXAMPLES:.ml=.byte.out) $(ALL_EXAMPLES:.ml=.opt.out)
	-@rm -f $(ALL_EXAMPLES:.ml=.sundials.out)
	-@rm -f $(ALL_EXAMPLES:.ml=.annot)
	-@rm -f $(ALL_EXAMPLES:.ml=.byte.diff) $(ALL_EXAMPLES:.ml=.opt.diff)
	-@rm -f $(ALL_EXAMPLES:.ml=.self.diff)
	-@rm -f tests.log tests.self.log tests.byte.log tests.opt.log

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

