include ../../../config

pendulum.byte: showpendulum.cmo pendulum.ml
	$(OCAMLC) $(OCAMLFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ -dllpath ../../../ \
	    bigarray.cma unix.cma graphics.cma \
	    sundials_ida.cma showpendulum.cmo pendulum.ml

pendulum.opt: showpendulum.cmx pendulum.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ \
	    $(INCLUDES) -I ../../../ \
	    bigarray.cmxa unix.cmxa graphics.cmxa \
	    sundials_ida.cmxa showpendulum.cmx pendulum.ml

showpendulum.cmo: showpendulum.ml
showpendulum.cmx: showpendulum.ml

clean:
	-rm -f showpendulum.cmi showpendulum.cmo showpendulum.cmx showpendulum.o
	-rm -f pendulum.cmi pendulum.cmo pendulum.cmx pendulum.o

realclean:cleanall
cleanall: clean
	-rm -f pendulum.byte pendulum.opt

