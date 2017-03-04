include ../../../config

SRCROOT = ../../../src

SUNDIALS_CMA = $(SRCROOT)/sundials.cma
SUNDIALS_CMXA = $(SRCROOT)/sundials.cmxa

EXAMPLES = arraydense.byte arraydense.opt 		\
	   dense.byte dense.opt				\
	   sparse.byte sparse.opt			\
	   arrayband.byte arrayband.opt 		\
	   band.byte band.opt				\
	   noband.byte noband.opt			\
	   spils_nvector1.byte spils_nvector1.opt	\
	   spils_serial1.byte spils_serial1.opt
CEXAMPLES = arraydense.exe dense.exe arrayband.exe band.exe sparse.exe

CFLAGS += -Wall

all: byte opt exe diff

byte: $(EXAMPLES)
opt: $(EXAMPLES:.byte=.opt)
exe: $(CEXAMPLES)
log: $(EXAMPLES:.byte=.byte.log) \
     $(EXAMPLES:.opt=.opt.log)   \
     $(CEXAMPLES:.exe=.exe.log)
diff: $(CEXAMPLES:.exe=.byte.diff) $(CEXAMPLES:.exe=.opt.diff)
	@echo ======================================================================
	@echo 'Summary (each should be 0):'
	@wc -l $^

arraydense.byte.diff: arraydense.byte.log
arraydense.opt.diff: arraydense.opt.log
dense.byte.diff: dense.byte.log
dense.opt.diff: dense.opt.log
arrayband.byte.diff: arrayband.byte.log
arrayband.opt.diff: arrayband.opt.log
band.byte.diff: band.byte.log
band.opt.diff: band.opt.log
sparse.byte.diff: sparse.byte.log
sparse.opt.diff: sparse.opt.log

arrayband.byte: arrayband.ml $(SUNDIALS_CMA)
arrayband.opt: arrayband.ml $(SUNDIALS_CMXA)

band.byte: band.ml $(SUNDIALS_CMA)
band.opt: band.ml $(SUNDIALS_CMXA)

noband.byte: noband.ml $(SUNDIALS_CMA)
noband.opt: noband.ml $(SUNDIALS_CMXA)

spils_nvector1.byte: spils_nvector1.ml $(SUNDIALS_CMA)
spils_nvector1.opt: spils_nvector1.ml $(SUNDIALS_CMXA)

spils_serial1.byte: spils_serial1.ml $(SUNDIALS_CMA)
spils_serial1.opt: spils_serial1.ml $(SUNDIALS_CMXA)

clean:
	-@rm -f $(EXAMPLES) $(CEXAMPLES)
	-@rm -f $(EXAMPLES:.byte=.cmo) $(EXAMPLES:.byte=.cmx)
	-@rm -f $(EXAMPLES:.byte=.cmt) $(EXAMPLES:.byte=.cmti)
	-@rm -f $(EXAMPLES:.byte=.o) $(EXAMPLES:.byte=.cmi)
	-@rm -f $(EXAMPLES:.byte=.annot)
	-@rm -f $(EXAMPLES:.byte=.byte.log) $(CEXAMPLES:.exe=.byte.diff)
	-@rm -f $(EXAMPLES:.byte=.opt.log) $(CEXAMPLES:.exe=.opt.diff)
	-@rm -f $(CEXAMPLES:.exe=.exe.log)

distclean: clean
	-@rm -f $(EXAMPLES) $(EXAMPLES:.byte=.opt) $(CEXAMPLES)

# #

.SUFFIXES : .ml .byte .opt .c .exe .byte.log .opt.log .exe.log \
    	    .byte.diff .opt.diff

.c.exe:
	cc -o $@ $(CFLAGS) $(CVODE_CFLAGS) $< $(LIB_PATH) $(CVODE_LDFLAGS)

.ml.byte:
	$(OCAMLC) $(OCAMLFLAGS) -o $@ $(INCLUDES) \
	    -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) -dllpath $(SRCROOT) \
	    bigarray.cma unix.cma sundials.cma $<

.ml.opt:
	$(OCAMLOPT) $(OCAMLOPTFLAGS) -o $@ $(INCLUDES) \
	    -I $(SRCROOT) $(SUBDIRS:%=-I $(SRCROOT)/%) \
	    bigarray.cmxa unix.cmxa sundials.cmxa $<

.byte.byte.log:
	./$< > $@
.opt.opt.log:
	./$< > $@
.exe.exe.log:
	./$< > $@

.exe.log.byte.diff:
	diff -u $< $(<:.exe.log=.byte.log) > $@ || true

.exe.log.opt.diff:
	diff -u $< $(<:.exe.log=.opt.log) > $@ || true

