include ../../../config

# Invoke `make' with USELIB=sundials to run the tests with the
# sensitivity-agnostic subset of IDAS.  Note that memory usage
# statistics will differ from the versions in sundials/C, unless those
# are recompiled to link against IDAS.
USELIB ?= sundials_nosensi

MPIRUN = mpirun
EXAMPLES = idaFoodWeb_kry_p.ml idaFoodWeb_kry_bbd_p.ml \
	   idaHeat2D_kry_p.ml idaHeat2D_kry_bbd_p.ml

DIVIDER = "----------------------------------------------------------------------"

##

all: byte opt
byte: $(EXAMPLES:.ml=.byte)
opt: $(EXAMPLES:.ml=.opt)

genlog = \
    @(for f in $(1); do  			\
	echo $(DIVIDER);                 	\
	echo "--ida/parallel/$$f";         	\
	cat $$f;                         	\
     done;                               	\
     echo "Summary (each should be 0):"; 	\
     for f in $(1); do  			\
	(cd ../..; wc -l ida/parallel/$$f);	\
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

idaFoodWeb_kry_p.byte.out: idaFoodWeb_kry_p.byte
	$(MPIRUN) -np 4 $< > $@
idaFoodWeb_kry_p.opt.out: idaFoodWeb_kry_p.opt
	$(MPIRUN) -np 4 $< > $@
idaFoodWeb_kry_p.sundials.out: $(EXAMPLESROOT)/ida/parallel/idaFoodWeb_kry_p
	$(MPIRUN) -np 4 $< > $@

idaFoodWeb_kry_bbd_p.byte.out: idaFoodWeb_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
idaFoodWeb_kry_bbd_p.opt.out: idaFoodWeb_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
idaFoodWeb_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/ida/parallel/idaFoodWeb_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

idaHeat2D_kry_bbd_p.byte.out: idaHeat2D_kry_bbd_p.byte
	$(MPIRUN) -np 4 $< > $@
idaHeat2D_kry_bbd_p.opt.out: idaHeat2D_kry_bbd_p.opt
	$(MPIRUN) -np 4 $< > $@
idaHeat2D_kry_bbd_p.sundials.out: $(EXAMPLESROOT)/ida/parallel/idaHeat2D_kry_bbd_p
	$(MPIRUN) -np 4 $< > $@

idaHeat2D_kry_p.byte.out: idaHeat2D_kry_p.byte
	$(MPIRUN) -np 4 $< > $@
idaHeat2D_kry_p.opt.out: idaHeat2D_kry_p.opt
	$(MPIRUN) -np 4 $< > $@
idaHeat2D_kry_p.sundials.out: $(EXAMPLESROOT)/ida/parallel/idaHeat2D_kry_p
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
	    $(USELIB).cma sundials_mpi.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ bigarray.cmxa unix.cmxa mpi.cmxa \
	    $(USELIB).cmxa sundials_mpi.cmxa $<

%.byte.out: %.byte
	./$< > $@

%.opt.out: %.opt
	./$< > $@

%.sundials.out : $(EXAMPLESROOT)/ida/parallel/%
	$< > $@

