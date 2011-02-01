
OCAMLFLAGS =
OCAMLOPTFLAGS =

CC = gcc
CFLAGS = # -DRESTRICT_INTERNAL_PRECISION"

AR = ar

OCAMLC = $(OCAML_BIN)ocamlc.opt
OCAMLFLAGS =
OCAML_INCLUDE=`${OCAMLC} -where`
OCAML_LIBLINK=

OCAMLOPT = $(OCAML_BIN)ocamlopt.opt
OCAMLOPTFLAGS =

OCAMLDEP = $(OCAML_BIN)ocamldep

GNUPLOT = gnuplot

LIB = /usr/local/lib

# If sundials is configured with --with-blas or --with-lapack
# then the extra library dependency must also be included below
LAPACK_LIB = -cclib -lSimTKlapack

UNIX = unix.cma
UNIXX = unix.cmxa

DIRECTORIES = 

INCLUDES = $(DIRECTORIES:%=-I %)

OBJ_OPT = $(OBJ:.cmo=.cmx)

# Common rules
.SUFFIXES : .mli .ml .cmi .cmo .cmx .lci .lsi .c .o .cma .cmxa

.ml.cmo:
	$(OCAMLC) $(OCAMLFLAGS) -c $(INCLUDES) $<

.mli.cmi:
	$(OCAMLC) $(OCAMLFLAGS) -c $(INCLUDES) $<

.ml.cmx:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -c $(INCLUDES) $<

.cmo.cma:
	$(OCAMLC) -a $(OCAMLFLAGS) \
	    -o $@ -custom $< $(OCAML_LIBLINK)

.cmx.cmxa:
	$(OCAMLOPT) -a $(OCAMLOPTFLAGS) \
	    -o $@ $< $(OCAML_LIBLINK)

.c.o:
	$(CC) -I $(OCAML_INCLUDE) $(CFLAGS) -o $@ -c $<

