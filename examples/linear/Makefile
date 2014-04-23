include ../../config

EXAMPLES = arraydense.byte arraydense \
	   dense.byte dense \
	   arrayband.byte arrayband \
	   band.byte band \
	   noband.byte

CVODE_LIBS := $(shell sundials-config -t s -l c -m cvode -s libs) -llapack
CFLAGS += -Wall

all: $(EXAMPLES)
opt: $(EXAMPLES:.byte=.opt)

arraydense.byte: arraydense.ml
arraydense: arraydense.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)
dense.byte: dense.ml
dense: dense.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)

arrayband.byte: arrayband.ml
arrayband: arrayband.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)
band.byte: band.ml
band: band.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)

noband.byte: noband.ml

clean:
	-rm -f $(EXAMPLES:.byte=.cmo) $(EXAMPLES:.byte=.cmx)
	-rm -f $(EXAMPLES:.byte=.o) $(EXAMPLES:.byte=.cmi)

realclean: cleanall
cleanall: clean
	-rm -f $(EXAMPLES) $(EXAMPLES:.byte=.opt)

# #

.SUFFIXES : .ml .byte .opt

.ml.byte:
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    $(INCLUDES) -I ../../ -dllpath ../../ \
	    bigarray.cma unix.cma sundials.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../ \
	    bigarray.cmxa unix.cmxa sundials.cmxa $<

