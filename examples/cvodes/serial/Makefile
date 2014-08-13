include ../../../config

EXAMPLES = cvsAdvDiff_ASAi_bnd.ml \
	   cvsAdvDiff_FSA_non.ml  \
	   cvsDiurnal_FSA_kry.ml  \
	   cvsFoodWeb_ASAi_kry.ml \
	   cvsFoodWeb_ASAp_kry.ml \
	   cvsHessian_ASA_FSA.ml  \
	   cvsRoberts_ASAi_dns.ml \
	   cvsRoberts_FSA_dns.ml

DIVIDER = "----------------------------------------------------------------------"

##

all: byte opt
byte: $(EXAMPLES:.ml=.byte)
opt: $(EXAMPLES:.ml=.opt)

genlog = \
    @(for f in $(1); do  			\
	echo $(DIVIDER);                 	\
	echo "--cvodes/serial/$$f";      	\
	cat $$f;                         	\
     done;                               	\
     echo "Summary (each should be 0):"; 	\
     for f in $(1); do  			\
	(cd ../..; wc -l cvodes/serial/$$f);	\
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

cvsAdvDiff_FSA_non.sundials.out:$(EXAMPLESROOT)/cvodes/serial/cvsAdvDiff_FSA_non
	$< -sensi sim t > $@

cvsDiurnal_FSA_kry.sundials.out:$(EXAMPLESROOT)/cvodes/serial/cvsDiurnal_FSA_kry
	$< -sensi sim t > $@

cvsRoberts_FSA_dns.sundials.out:$(EXAMPLESROOT)/cvodes/serial/cvsRoberts_FSA_dns
	$< -sensi sim t > $@

##

distclean: clean
clean:
	-@rm -f $(EXAMPLES:.ml=.cmo)   $(EXAMPLES:.ml=.cmx)
	-@rm -f $(EXAMPLES:.ml=.o)     $(EXAMPLES:.ml=.cmi)
	-@rm -f $(EXAMPLES:.ml=.c.log) $(EXAMPLES:.ml=.ml.log)
	-@rm -f $(EXAMPLES:.ml=.byte) $(EXAMPLES:.ml=.opt)
	-@rm -f $(EXAMPLES:.ml=.byte.out) $(EXAMPLES:.ml=.opt.out)
	-@rm -f $(EXAMPLES:.ml=.sundials.out)
	-@rm -f $(EXAMPLES:.ml=.annot)
	-@rm -f $(EXAMPLES:.ml=.byte.diff) $(EXAMPLES:.ml=.opt.diff)
	-@rm -f $(EXAMPLES:.ml=.self.diff)
	-@rm -f tests.log tests.self.log
	-@rm -f tests.byte.log tests.opt.log

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

