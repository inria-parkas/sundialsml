SRCROOT=../../..
SUBDIR=cvodes/serial

EXAMPLES = cvsAdvDiff_ASAi_bnd.ml \
	   cvsAdvDiff_FSA_non.ml  \
	   cvsDiurnal_FSA_kry.ml  \
	   cvsFoodWeb_ASAi_kry.ml \
	   cvsFoodWeb_ASAp_kry.ml \
	   cvsHessian_ASA_FSA.ml  \
	   cvsRoberts_ASAi_dns.ml \
	   cvsRoberts_FSA_dns.ml

include ../../examples.mk

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
