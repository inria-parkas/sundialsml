SRCROOT = ../../..
SUBDIR = idas/serial

EXAMPLES = idasRoberts_FSA_dns.ml idasAkzoNob_ASAi_dns.ml idasAkzoNob_dns.ml \
           idasHessian_ASA_FSA.ml idasRoberts_ASAi_dns.ml idasSlCrank_FSA_dns.ml

include ../../examples.mk

# Tests with command-line arguments (make same choice as Sundials distrib)

idasRoberts_FSA_dns.byte.out: idasRoberts_FSA_dns.byte
	./$< -sensi sim t > $@

idasRoberts_FSA_dns.opt.out: idasRoberts_FSA_dns.opt
	./$< -sensi sim t > $@

idasRoberts_FSA_dns.sundials.out: $(EXAMPLESROOT)/idas/serial/idasRoberts_FSA_dns
	$< -sensi sim t > $@
