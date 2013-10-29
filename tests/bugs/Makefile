include ../../config

TESTS=cvode_stop_time.byte ida_stop_time.byte

SUNDIALS=-I ../.. unix.cma bigarray.cma sundials.cma

all: $(TESTS)

%.byte: %.ml
	$(OCAMLC) $(OCAMLFLAGS) -o $@ $(SUNDIALS) $<

clean:
	-rm -f $(TESTS:.byte=.cmo)
	-rm -f $(TESTS:.byte=.cmi)
	-rm -f $(TESTS)
