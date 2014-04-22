include ../../config

EXAMPLES = arraydense.byte arraydense dense

CVODE_LIBS := $(shell sundials-config -t s -l c -m cvode -s libs) -llapack
CFLAGS += -Wall

all: $(EXAMPLES)
opt: $(EXAMPLES:.byte=.opt)

arraydense.byte: arraydense.ml
arraydense: arraydense.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)
dense: dense.c
	cc -o $@ $(CFLAGS) $< $(CVODE_LIBS)

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

