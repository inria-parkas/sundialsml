include ../../../config

SUNDIALS_CMA = ../../../sundials.cma
SUNDIALS_CMXA = ../../../sundials.cmxa

EXAMPLES = arraydense.byte arraydense.opt arraydense \
	   dense.byte dense.opt dense \
	   arrayband.byte arrayband.opt arrayband \
	   band.byte band.opt band \
	   noband.byte noband.opt \
	   spils_nvector1.byte spils_nvector1.opt \
	   spils_serial1.byte spils_serial1.opt

CVODE_LIBS := $(shell sundials-config -t s -l c -m cvode -s libs) -llapack
CFLAGS += -Wall

all: byte opt

byte: $(EXAMPLES)
opt: $(EXAMPLES:.byte=.opt)

arraydense.byte: arraydense.ml
arraydense.opt: arraydense.ml
arraydense: arraydense.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)
dense.byte: dense.ml
dense.opt: dense.ml
dense: dense.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)

arrayband.byte: arrayband.ml $(SUNDIALS_CMA)
arrayband.opt: arrayband.ml $(SUNDIALS_CMXA)
arrayband: arrayband.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)
band.byte: band.ml $(SUNDIALS_CMA)
band.opt: band.ml $(SUNDIALS_CMXA)
band: band.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)

noband.byte: noband.ml $(SUNDIALS_CMA)
noband.opt: noband.ml $(SUNDIALS_CMXA)

spils_nvector1.byte: spils_nvector1.ml $(SUNDIALS_CMA)
spils_nvector1.opt: spils_nvector1.ml $(SUNDIALS_CMXA)

spils_serial1.byte: spils_serial1.ml $(SUNDIALS_CMA)
spils_serial1.opt: spils_serial1.ml $(SUNDIALS_CMXA)

clean:
	-@rm -f $(EXAMPLES:.byte=.cmo) $(EXAMPLES:.byte=.cmx)
	-@rm -f $(EXAMPLES:.byte=.o) $(EXAMPLES:.byte=.cmi)
	-@rm -f $(EXAMPLES:.byte=.annot)

distclean: clean
	-@rm -f $(EXAMPLES) $(EXAMPLES:.byte=.opt)

# #

.SUFFIXES : .ml .byte .opt

.ml.byte:
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ -dllpath ../../../ \
	    bigarray.cma unix.cma sundials.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ \
	    bigarray.cmxa unix.cmxa sundials.cmxa $<

