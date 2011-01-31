include ../../Makefile.inc

ball.byte: showball.cmo ball.ml
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    $(INCLUDES) -I ../../ \
	    bigarray.cma unix.cma graphics.cma \
	    ../../cvode.cma showball.cmo ball.ml

ball.opt: showball.cmx ball.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../ \
	    bigarray.cmxa unix.cmxa graphics.cmxa \
	    ../../cvode.cmxa showball.cmx ball.ml

showball.cmo: showball.cmi showball.ml
showball.cmx: showball.cmi showball.ml
showball.cmi: showball.mli

clean:
	-rm -f showball.cmi showball.cmo showball.cmx showball.o
	-rm -f ball.cmi ball.cmo ball.cmx ball.o

cleanall: clean
	-rm -f ball.byte ball.opt

