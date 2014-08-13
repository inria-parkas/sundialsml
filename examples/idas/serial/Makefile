include ../../../config

EXAMPLES = idasRoberts_FSA_dns.ml idasAkzoNob_ASAi_dns.ml idasAkzoNob_dns.ml \
           idasHessian_ASA_FSA.ml idasRoberts_ASAi_dns.ml idasSlCrank_FSA_dns.ml

DIVIDER = "----------------------------------------------------------------------"

##

all: byte opt
byte: $(EXAMPLES:.ml=.byte)
opt: $(EXAMPLES:.ml=.opt)

genlog = \
    @(for f in $(1); do  			\
	echo $(DIVIDER);                 	\
	echo "--idas/serial/$$f";        	\
	cat $$f;                         	\
     done;                               	\
     echo "Summary (each should be 0):"; 	\
     for f in $(1); do  			\
	(cd ../..; wc -l idas/serial/$$f);	\
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

idasRoberts_FSA_dns.opt.out: idasRoberts_FSA_dns.opt
	./$< -sensi sim t > $@

idasRoberts_FSA_dns.sundials.out:$(EXAMPLESROOT)/idas/serial/idasRoberts_FSA_dns
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

%.sundials.out : $(EXAMPLESROOT)/idas/serial/%
	$< > $@

