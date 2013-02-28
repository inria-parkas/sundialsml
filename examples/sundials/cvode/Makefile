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

all: $(EXAMPLES:.ml=.byte)
opt: $(EXAMPLES:.ml=.opt)

lapack: $(LAPACK_EXAMPLES:.ml=.byte)
lapack.opt: $(LAPACK_EXAMPLES:.ml=.opt)

tests.log: $(EXAMPLES:.ml=.diff)
	@(for f in $(EXAMPLES:.ml=.diff); do  \
	    echo $(DIVIDER);                  \
	    echo "--$$f";                    \
	    cat $$f;                         \
	 done;                               \
	 echo "Summary (each should be 0):"; \
	 for f in $(EXAMPLES:.ml=.diff); do  \
	    wc -l $$f;                       \
	 done;                               \
	) > $@

lapack-tests.log: $(LAPACK_EXAMPLES:.ml=.diff)
	@(for f in $(LAPACK_EXAMPLES:.ml=.diff); do  \
	    echo $(DIVIDER);                  \
	    echo "--$$f";                    \
	    cat $$f;                         \
	 done;                               \
	 echo "Summary (each should be 0):"; \
	 for f in $(LAPACK_EXAMPLES:.ml=.diff); do  \
	    wc -l $$f;                       \
	 done;                               \
	) > $@

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
	-@rm -f $(ALL_EXAMPLES:.ml=.diff) tests.log lapack-tests.log

# #

.SUFFIXES : .ml .byte .opt .diff

.ml.byte:
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ -dllpath ../../../ bigarray.cma unix.cma \
	    sundials_cvode.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ bigarray.cmxa unix.cmxa \
	    sundials_cvode.cmxa $<

.byte.diff:
	cat $(EXAMPLESROOT)/cvode/serial/$(<:.byte=.out) > $(@:.diff=.c.log)
	./$< > $(@:.diff=.ml.log)
	diff -u $(@:.diff=.c.log) $(@:.diff=.ml.log) > $@ || true

