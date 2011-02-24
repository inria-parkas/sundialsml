
OCAMLFLAGS =
OCAMLOPTFLAGS =

CC = gcc
CFLAGS = # -DRESTRICT_INTERNAL_PRECISION"

AR = ar

OCAMLC = $(OCAML_BIN)ocamlc.opt
OCAMLFLAGS =
OCAML_INCLUDE := $(shell $(OCAMLC) -where)
OCAML_LIBLINK=

PKGDIR = $(OCAML_INCLUDE)/sundials

OCAMLOPT = $(OCAML_BIN)ocamlopt.opt
OCAMLOPTFLAGS =

OCAMLDEP = $(OCAML_BIN)ocamldep

OCAMLMKLIB = $(OCAML_BIN)ocamlmklib -ocamlc "$(OCAMLC)" -ocamlopt "$(OCAMLOPT)"
OCAMLMKLIBFLAGS = -verbose

GNUPLOT = gnuplot

LIB = /usr/local/lib

uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
ifneq (,$(findstring MINGW,$(uname_S)))
XA=.lib
XO=.obj
XS=.dll
else
XA=.a
XO=.o
XS=.so
endif

CP=cp
RM=rm
CAT=cat
ECHO=echo
MKDIR=mkdir -p
RMDIR=rmdir

#EXAMPLES_ROOT = /usr/local/examples
EXAMPLES_ROOT = /home/tbourke/software/sundials-2.4.0/examples

# If sundials is configured with --with-blas or --with-lapack
# then the extra library dependency must also be included below
LAPACK_LIB = -lSimTKlapack

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

