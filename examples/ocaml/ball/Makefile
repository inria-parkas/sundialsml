include ../../../config

SRCROOT = ../../../src

all: ball.byte ball.opt

ball.byte: showball.cmi showball.cmo ball.ml
	$(OCAMLC) $(OCAMLFLAGS) -o $@ $(INCLUDES) \
	    -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) -dllpath $(SRCROOT) \
	    unix.cma bigarray.cma graphics.cma \
	    sundials.cma showball.cmo ball.ml

ball.opt: showball.cmi showball.cmx ball.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ $(INCLUDES) \
	    -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) \
	    unix.cmxa bigarray.cmxa graphics.cmxa \
	    sundials.cmxa showball.cmx ball.ml

showball.cmi : showball.mli
	$(OCAMLC) $(OCAMLFLAGS) $^

showball.cmo : showball.ml
	$(OCAMLC) $(OCAMLFLAGS) unix.cma graphics.cma -c $^

showball.cmx : showball.ml
	$(OCAMLOPT) $(OCAMLOPTFLAGS) unix.cmxa graphics.cmxa -c $^

clean:
	-@rm -f showball.cmi showball.cmo showball.cmx showball.o
	-@rm -f ball.cmi ball.cmo ball.cmx ball.o
	-@rm -f showball.annot ball.annot

distclean: clean
	-@rm -f ball.byte ball.opt

