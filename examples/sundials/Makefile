include ../../Makefile.inc

EXAMPLES = cvAdvDiff_bnd.ml \
	   cvRoberts_dns.ml \
	   cvRoberts_dns_uw.ml \
	   cvDirectDemo_ls.ml \
	   cvDiurnal_kry.ml \
	   cvDiurnal_kry_bp.ml

LAPACK_EXAMPLES = cvAdvDiff_bndL.ml \
	   	  cvRoberts_dnsL.ml

ALL_EXAMPLES = $(EXAMPLES) $(LAPACK_EXAMPLES)

all: $(EXAMPLES:.ml=.byte)
opt: $(EXAMPLES:.ml=.opt)

lapack: $(LAPACK_EXAMPLES:.ml=.byte)
lapack.opt: $(LAPACK_EXAMPLES:.ml=.opt)

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

clean:
	-rm -f $(ALL_EXAMPLES:.ml=.cmo) $(ALL_EXAMPLES:.ml=.cmx)
	-rm -f $(ALL_EXAMPLES:.ml=.o) $(ALL_EXAMPLES:.ml=.cmi)

cleanall: clean
	-rm -f $(ALL_EXAMPLES:.ml=.byte) $(ALL_EXAMPLES:.ml=.opt)

# #

.SUFFIXES : .ml .byte .opt

.ml.byte:
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    $(INCLUDES) -I ../../ bigarray.cma unix.cma ../../cvode.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../ bigarray.cmxa unix.cmxa ../../cvode.cmxa $<

