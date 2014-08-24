include ../../../config

# The following examples are omitted because they are identical to their
# counterparts in ida (without `s'):
# idasHeat2D_kry_p.c
# idasHeat2D_kry_bbd_p.c
# idasFoodWeb_kry_p.c
# idasFoodWeb_kry_bbd_p.c

MPIRUN = mpirun
EXAMPLES = idasHeat2D_FSA_kry_bbd_p.ml		\
	   idasBruss_kry_bbd_p.ml		\
	   idasBruss_FSA_kry_bbd_p.ml		\
	   idasBruss_ASAp_kry_bbd_p.ml

DIVIDER = "----------------------------------------------------------------------"

##

all: byte opt
byte: $(EXAMPLES:.ml=.byte)
opt: $(EXAMPLES:.ml=.opt)

genlog = \
    @(for f in $(1); do  			\
	echo $(DIVIDER);                 	\
	echo "--idas/parallel/$$f";    	\
	cat $$f;                         	\
     done;                               	\
     echo "Summary (each should be 0):"; 	\
     for f in $(1); do  			\
	(cd ../..; wc -l idas/parallel/$$f);	\
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
	@grep -q '^[^0]' $@

tests.opt.log: $(EXAMPLES:.ml=.opt.out) $(EXAMPLES:.ml=.sundials.out)
	@for f in $(EXAMPLES:.ml=); do \
	    diff -u $$f.sundials.out $$f.opt.out > $$f.opt.diff || true; \
	done
	$(call genlog, $(EXAMPLES:.ml=.opt.diff), $@)
	cat $@
	@grep -q '^[^0]' $@

tests.self.log: $(EXAMPLES:.ml=.byte.out) $(EXAMPLES:.ml=.opt.out)
	@for f in $(EXAMPLES:.ml=); do \
	    diff -u $$f.byte.out $$f.opt.out > $$f.self.diff || true; \
	done
	$(call genlog, $(EXAMPLES:.ml=.self.diff), $@)
	cat $@
	@grep -q '^[^0]' $@

##

# This test cannot be run with -sensi sim.
idasHeat2D_FSA_kry_bbd_p.byte.out: idasHeat2D_FSA_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< -sensi stg t > $@
idasHeat2D_FSA_kry_bbd_p.opt.out: idasHeat2D_FSA_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< -sensi stg t > $@
idasHeat2D_FSA_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/idas/parallel/idasHeat2D_FSA_kry_bbd_p
	$(MPIRUN) -np 4 $< -sensi stg t > $@

idasBruss_kry_bbd_p.byte.out: idasBruss_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
idasBruss_kry_bbd_p.opt.out: idasBruss_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
idasBruss_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/idas/parallel/idasBruss_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

idasBruss_FSA_kry_bbd_p.byte.out: idasBruss_FSA_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
idasBruss_FSA_kry_bbd_p.opt.out: idasBruss_FSA_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
idasBruss_FSA_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/idas/parallel/idasBruss_FSA_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

idasBruss_ASAp_kry_bbd_p.byte.out: idasBruss_ASAp_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
idasBruss_ASAp_kry_bbd_p.opt.out: idasBruss_ASAp_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
idasBruss_ASAp_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/idas/parallel/idasBruss_ASAp_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

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
	    $(INCLUDES) -I ../../../ -dllpath ../../../ \
	    bigarray.cma unix.cma mpi.cma \
	    sundials.cma sundials_mpi.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ bigarray.cmxa unix.cmxa mpi.cmxa \
	    sundials.cmxa sundials_mpi.cmxa $<

%.byte.out: %.byte
	./$< > $@

%.opt.out: %.opt
	./$< > $@

%.sundials.out : $(EXAMPLESROOT)/ida/parallel/%
	$< > $@

