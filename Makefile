include Makefile.inc

MLOBJ = nvector.cmo \
	cvode.cmo \
	nvector_array.cmo

COBJ =	ml_cvode.o \
	ml_cvode_bp.o \
	ml_cvode_ba.o \
	ml_cvode_nvec.o \
	ml_nvector.o \

# TODO: Remove solvelucy and examples that depend on it
#       (just keep the examples converted from Sundials) 

all: cvode.cma cvode.cmxa nvector.cma nvector.cmxa \
    	nvector_array.cma nvector_array.cmxa

cvode.cma: OCAML_LIBLINK := -cclib -lsundials_cvode -cclib \
    			    -lsundials_nvecserial \
    			    $(LAPACK_LIB) -cclib -lmlcvode
cvode.cma: libmlcvode.a

cvode.cmxa: OCAML_LIBLINK := -cclib -lsundials_cvode -cclib \
    			     -lsundials_nvecserial \
    			     $(LAPACK_LIB) -cclib -lmlcvode
cvode.cmxa: libmlcvode.a

libmlcvode.a: ml_cvode.o ml_cvode_bp.o \
    		ml_cvode_nvec.o ml_cvode_ba.o ml_nvector.o
	$(AR) rc libmlcvode.a ml_cvode.o ml_cvode_bp.o \
		 ml_cvode_nvec.o ml_cvode_ba.o ml_nvector.o

nvector.cma: OCAML_LIBLINK := -cclib -lsundials_cvode -cclib -lmlcvode
nvector.cma: nvector.cmi libmlcvode.a

nvector_array.cma: OCAML_LIBLINK := -cclib -lsundials_cvode -cclib -lmlcvode
nvector_array.cma: nvector_array.cmi nvector.cmi libmlcvode.a

nvector.cmxa: OCAML_LIBLINK := -cclib -lsundials_cvode -cclib -lmlcvode
nvector.cmxa: nvector.cmi libmlcvode.a

nvector_array.cmxa: OCAML_LIBLINK := -cclib -lsundials_cvode -cclib -lmlcvode
nvector_array.cmxa: nvector_array.cmi nvector.cmi libmlcvode.a

ml_cvode.o: ml_cvode.c
ml_cvode_bp.o: ml_cvode_bp.c
ml_cvode_nvec.o: ml_cvode_nvec.c
ml_nvector.o: ml_nvector.c

ml_cvode_ba.o: ml_cvode_nvec.c
	$(CC) -I $(OCAML_INCLUDE) $(CFLAGS) \
	    -DML_CVODE_BIGARRAYS -o $@ -c $<

# ##

depend:
	$(OCAMLDEP) $(INCLUDES) *.mli *.ml > .depend

# TODO
clean:
	-rm -f $(MLOBJ) $(MLOBJ:.cmo=.cmx)
	-rm -f $(COBJ)
	-rm -f cvode.o nvector.o nvector_array.o

# TODO
cleanall: clean
	-rm -f libmlcvode.a cvode.a
	-rm -f nvector.a nvector_array.a
	-rm -f cvode.cma cvode.cmxa
	-rm -f nvector.cma nvector.cmxa
	-rm -f nvector_array.cma nvector_array.cmxa
	-rm -f cvode.cmi nvector.cmi nvector_array.cmi

-include .depend
