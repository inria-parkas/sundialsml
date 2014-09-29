include ../../../config

EXAMPLES = arraydense.byte arraydense.opt arraydense \
	   dense.byte dense.opt dense2.byte dense2.opt dense \
	   arrayband.byte arrayband.opt arrayband \
	   band.byte band.opt band2.byte band2.opt band \
	   noband.byte noband.opt \
	   spils_nvector1.byte spils_nvector1.opt \
	   spils_serial1.byte spils_serial1.opt

CVODE_LIBS := $(shell sundials-config -t s -l c -m cvode -s libs) -llapack
CFLAGS += -Wall

all: $(EXAMPLES)
opt: $(EXAMPLES:.byte=.opt)

arraydense.byte: arraydense.ml
arraydense.opt: arraydense.ml
arraydense: arraydense.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)
dense.byte: dense.ml
dense.opt: dense.ml
dense2.byte: dense2.ml
dense2.opt: dense2.ml
dense: dense.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)

arrayband.byte: arrayband.ml
arrayband.opt: arrayband.ml
arrayband: arrayband.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)
band.byte: band.ml
band.opt: band.ml
band2.byte: band2.ml
band2.opt: band2.ml
band: band.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)

noband.byte: noband.ml
noband.opt: noband.ml

spils_nvector1.byte: spils_nvector1.ml
spils_nvector1.opt: spils_nvector1.ml

spils_serial1.byte: spils_serial1.ml
spils_serial1.opt: spils_serial1.ml

clean:
	-@rm -f $(EXAMPLES:.byte=.cmo) $(EXAMPLES:.byte=.cmx)
	-@rm -f $(EXAMPLES:.byte=.o) $(EXAMPLES:.byte=.cmi)

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

