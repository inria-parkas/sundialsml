include ../../../config

EXAMPLES = cvAdvDiff_bnd.ml \
	   cvRoberts_dns.ml \
	   cvRoberts_dns_uw.ml \
	   cvDirectDemo_ls.ml \
	   cvDiurnal_kry.ml \
	   cvDiurnal_kry_bp.ml \
	   cvKrylovDemo_ls.ml \
	   cvKrylovDemo_prec.ml

LAPACK_EXAMPLES = cvAdvDiff_bndL.ml \
		  cvRoberts_dnsL.ml

ALL_EXAMPLES = $(EXAMPLES) $(LAPACK_EXAMPLES)

DIVIDER = "----------------------------------------------------------------------"

##

all: byte opt
byte: $(EXAMPLES:.ml=.byte)
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

cvAdvDiff_bnd.byte: cvAdvDiff_bnd.ml
cvAdvDiff_bnd.opt:  cvAdvDiff_bnd.ml

cvAdvDiff_bndL.byte: cvAdvDiff_bndL.ml
cvAdvDiff_bndL.opt:  cvAdvDiff_bndL.ml

cvRoberts_dns.byte: cvRoberts_dns.ml
cvRoberts_dns.opt:  cvRoberts_dns.ml

cvRoberts_dnsL.byte: cvRoberts_dnsL.ml
cvRoberts_dnsL.opt:  cvRoberts_dnsL.ml

cvRoberts_dns_uw.byte: cvRoberts_dns_uw.ml
cvRoberts_dns_uw.opt:  cvRoberts_dns_uw.ml

cvDirectDemo_ls.byte: cvDirectDemo_ls.ml
cvDirectDemo_ls.opt:  cvDirectDemo_ls.ml

cvDiurnal_kry.byte: cvDiurnal_kry.ml
cvDiurnal_kry.opt:  cvDiurnal_kry.ml

cvDiurnal_kry_bp.byte: cvDiurnal_kry_bp.ml
cvDiurnal_kry_bp.opt:  cvDiurnal_kry_bp.ml

cvKrylovDemo_ls.byte: cvKrylovDemo_ls.ml
cvKrylovDemo_ls.opt:  cvKrylovDemo_ls.ml

cvKrylovDemo_prec.byte: cvKrylovDemo_prec.ml
cvKrylovDemo_prec.opt:  cvKrylovDemo_prec.ml

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

%.sundials.out : $(EXAMPLESROOT)/cvode/serial/%
	$< > $@

