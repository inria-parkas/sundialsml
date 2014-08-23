include ../../../config

# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of CVODES.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against CVODES.
USELIB ?= sundials_nosensi

MPIRUN = mpirun

EXAMPLES = kinFoodWeb_kry_p.ml \
	   kinFoodWeb_kry_bbd_p.ml

ALL_EXAMPLES = $(EXAMPLES)

DIVIDER = "----------------------------------------------------------------------"

##

all: byte opt

byte: $(ALL_EXAMPLES:.ml=.byte)
opt: $(ALL_EXAMPLES:.ml=.opt)

genlog = \
    @(for f in $(1); do  			\
	echo $(DIVIDER);                 	\
	echo "--kinsol/parallel/$$f";    	\
	cat $$f;                         	\
     done;                               	\
     echo "Summary (each should be 0):"; 	\
     for f in $(1); do  			\
	(cd ../..; wc -l kinsol/parallel/$$f);	\
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

kinFoodWeb_kry_p.byte: kinFoodWeb_kry_p.ml
kinFoodWeb_kry_p.opt: kinFoodWeb_kry_p.ml
kinFoodWeb_kry_p.byte.out: kinFoodWeb_kry_p.byte
	$(MPIRUN) -np 4 $< > $@
kinFoodWeb_kry_p.opt.out: kinFoodWeb_kry_p.opt
	$(MPIRUN) -np 4 $< > $@
kinFoodWeb_kry_p.sundials.out: $(EXAMPLESROOT)/kinsol/parallel/kinFoodWeb_kry_p
	$(MPIRUN) -np 4 $< > $@

kinFoodWeb_kry_bbd_p.byte: kinFoodWeb_kry_bbd_p.ml
kinFoodWeb_kry_bbd_p.opt: kinFoodWeb_kry_bbd_p.ml
kinFoodWeb_kry_bbd_p.byte.out: kinFoodWeb_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
kinFoodWeb_kry_bbd_p.opt.out: kinFoodWeb_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
kinFoodWeb_kry_bbd_p.sundials.out: \
    $(EXAMPLESROOT)/kinsol/parallel/kinFoodWeb_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

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
	    $(INCLUDES) -I ../../../ -dllpath ../../../ \
	    bigarray.cma unix.cma mpi.cma \
	    $(USELIB).cma sundials_mpi.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ bigarray.cmxa unix.cmxa mpi.cmxa \
	    $(USELIB).cmxa sundials_mpi.cmxa $<

%.byte.out: %.byte
	./$< > $@

%.opt.out: %.opt
	./$< > $@

%.sundials.out : $(EXAMPLESROOT)/kinsol/parallel/%
	$< > $@

