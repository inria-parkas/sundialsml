include ../../../config

SRCROOT = ../../../src

all: sincos.byte sincos_nvec.byte sincos.opt sincos_nvec.opt

NVECTOR=

byte.plot: sincos.byte
	./$< | head -n 500 | grep -v '^R: ' > sincos.log
	./$< | head -n 500 | sed -ne 's/^R: //p' > sincos-roots.log
	gnuplot -persist sincos.gplot

opt.plot: sincos.opt
	./$< | head -n 500 | grep -v '^R: ' > sincos.log
	./$< | head -n 500 | sed -ne 's/^R: //p' > sincos-roots.log
	gnuplot -persist sincos.gplot

nvec-byte.plot: sincos_nvec.byte
	./$< | head -n 500 | grep -v '^R: ' > sincos.log
	./$< | head -n 500 | sed -ne 's/^R: //p' > sincos-roots.log
	gnuplot -persist sincos.gplot

nvec-opt.plot: sincos_nvec.opt
	./$< | head -n 500 | grep -v '^R: ' > sincos.log
	./$< | head -n 500 | sed -ne 's/^R: //p' > sincos-roots.log
	gnuplot -persist sincos.gplot

sincos.byte: sincos.ml

sincos_nvec.byte: sincos_nvec.ml

sincos.opt: sincos.ml

sincos_nvec.opt: sincos_nvec.ml

clean:
	-@rm -f sincos.cmi sincos.cmo sincos.cmx sincos.o
	-@rm -f sincos.cmt sincos.cmti
	-@rm -f sincos_nvec.cmi sincos_nvec.cmo sincos_nvec.cmx sincos_nvec.o
	-@rm -f sincos_nvec.cmt sincos_nvec.cmti
	-@rm -f sincos.log sincos-roots.log

cleanall distclean: clean
	-@rm -f sincos.byte sincos.opt
	-@rm -f sincos_nvec.opt sincos_nvec.byte

.SUFFIXES : .ml .byte .opt

.ml.byte:
	$(OCAMLC) $(OCAMLFLAGS) -o $@ $(INCLUDES) $(EXAMPLE_INCLUDES) \
	    -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) -dllpath $(SRCROOT) \
	    $(BIGARRAY_CMA) sundials.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ $(INCLUDES) $(EXAMPLE_INCLUDES) \
	    -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) $(LIB_PATH:%=-ccopt %) \
	    $(BIGARRAY_CMXA) sundials.cmxa $<

