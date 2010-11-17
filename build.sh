
# todo:
# * use a proper build system
# * also allow for native code
# * also allow for generation and dynamic linking of cvode.so

CC=cc
CFLAGS="" # -DRESTRICT_INTERNAL_PRECISION"
AR=ar
OCAMLC=ocamlc
OCAMLOPT=ocamlopt
OCAMLOPTFLAGS=""
LIB=/usr/local/lib
OCAML_INCLUDE=`${OCAMLC} -where`
GNUPLOT=gnuplot
MAX_LINES=500

# If sundials is configured with --with-blas or --with-lapack
# then the extra library dependency must also be included below
LAPACK_LIB= # "-cclib -lSimTKlapack"

BASIC_EXAMPLES="discontinuous sincos cchatter"
LUCYSOLVE_EXAMPLES="nontordu \
		    nontordu2 \
		    nontordu3 \
		    sincos_lucyf \
		    billiard1d \
		    cascade \
		    example3"
SUNDIALS_EXAMPLES="cvRoberts_dns cvAdvDiff_bnd"
SUNDIALS_NVECTOR_EXAMPLES="cvRoberts_dns_nvec"
PLOT="examples/billiard1d \
      examples/nontordu \
      examples/nontordu2 \
      examples/nontordu3 \
      examples/example3 " 

function plot_example
{
    if expr "$PLOT" : ".*/$1 .*" > /dev/null; then
	echo "* plotting: $f.ps"
	./$f | head -n ${MAX_LINES} > $f.log
	${GNUPLOT} -persist $f.gplot > $f.ps
    fi
}

case $1 in
clean)
    rm -f ml_cvode.o ml_cvode_bp.o libmlcvode.a
    rm -f ml_cvode_nvec.o ml_cvode_ba.o
    rm -f cvode.cmi cvode.cmo cvode.cma
    rm -f cvode.cmx cvode.cmxa
    rm -f nvector.cmo nvector.cmi ml_nvector.o
    rm -f nvector_array.cmo nvector_array.cmi
    rm -f solvelucy.cmi solvelucy.cmo

    rm -f solvelucy.cmx solvelucy.o

    rm -f examples/ball.cmi examples/ball.cmo examples/ball.cmx
    rm -f examples/ball.opt examples/ball.o examples/ball
    rm -f examples/showball.cmi examples/showball.cmo examples/showball.cma
    rm -f examples/showball.cmx examples/showball.o examples/showball.opt
    rm -f examples/showball.a
    rm -f examples/sincos.cmi examples/sincos.cmo
    rm -f examples/sincos_lucyf.cmi examples/sincos_lucyf.cmo
    rm -f examples/sincos examples/sincos_lucyf

    rm -f examples/nontordu examples/nontordu.cmo
    rm -f examples/nontordu.cmi

    for f in $PLOT; do
	rm -f $f.log
	rm -f $f.ps
    done

    for f in $BASIC_EXAMPLES $LUCYSOLVE_EXAMPLES; do
	rm -f examples/$f.cmo
	rm -f examples/$f.o
	rm -f examples/$f.cmx
	rm -f examples/$f.cmi
	rm -f examples/$f
	rm -f examples/$f.opt
    done

    for f in $SUNDIALS_EXAMPLES; do
	rm -f examples/sundials/$f.cmo
	rm -f examples/sundials/$f.o
	rm -f examples/sundials/$f.cmi
	rm -f examples/sundials/$f.cmx
	rm -f examples/sundials/$f
	rm -f examples/sundials/$f.opt
    done
    ;;

*)
    echo "* ml_cvode.c -> ml_cvode.o"
    ${CC} -I $OCAML_INCLUDE ${CFLAGS} -c ml_cvode.c || exit 1

    echo "* ml_cvode_bp.c -> ml_cvode_bp.o"
    ${CC} -I $OCAML_INCLUDE ${CFLAGS} -c ml_cvode_bp.c || exit 1

    echo "* ml_cvode_nvec.c -> ml_cvode_nvec.o"
    ${CC} -I $OCAML_INCLUDE ${CFLAGS} -c ml_cvode_nvec.c || exit 1

    echo "* ml_cvode_nvec.c -> ml_cvode_ba.o"
    ${CC} -I $OCAML_INCLUDE ${CFLAGS} -DML_CVODE_BIGARRAYS \
	-o ml_cvode_ba.o -c ml_cvode_nvec.c || exit 1

    echo "* ml_nvector.c -> ml_nvector.o"
    ${CC} -I $OCAML_INCLUDE ${CFLAGS} -o ml_nvector.o -c ml_nvector.c || exit 1

    echo "* ... -> libmlcvode.a"
    ${AR} rc libmlcvode.a \
	ml_cvode.o ml_cvode_bp.o ml_cvode_nvec.o ml_cvode_ba.o \
	ml_nvector.o || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* cvode.ml -> cvode.cmx"
	${OCAMLOPT} -c ${OCAMLOPTFLAGS} cvode.ml || exit 1
    fi

    echo "* ... -> cvode.cma"
    ${OCAMLC} -a -o cvode.cma -custom cvode.cmo \
	-cclib -lsundials_cvode \
	-cclib -lsundials_nvecserial \
	${LAPACK_LIB} \
	-cclib -lmlcvode || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* ... -> cvode.cmxa"
	${OCAMLOPT} -a -o cvode.cmxa ${OCAMLOPTFLAGS} \
	    cvode.cmx \
	    -cclib -lsundials_cvode \
	    -cclib -lsundials_nvecserial \
	    ${LAPACK_LIB} \
	    -cclib -lmlcvode || exit 1
    fi

    echo "* nvector.mli -> nvector.cmi"
    ${OCAMLC} nvector.mli || exit 1

    echo "* nvector.ml -> nvector.cmo"
    ${OCAMLC} -c nvector.ml || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* nvector.ml -> nvector.cmx"
	${OCAMLOPT} -c ${OCAMLOPTFLAGS} nvector.ml || exit 1
    fi

    echo "* ... -> nvector.cma"
    ${OCAMLC} -a -o nvector.cma -custom nvector.cmo \
	-cclib -lsundials_cvode \
	-cclib -lmlcvode || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* ... -> nvector.cmxa"
	${OCAMLOPT} -a -o nvector.cmxa ${OCAMLOPTFLAGS} \
	    nvector.cmx \
	    -cclib -lsundials_cvode \
	    -cclib -lmlcvode || exit 1
    fi

    echo "* nvector_array.mli -> nvector_array.cmi"
    ${OCAMLC} nvector_array.mli || exit 1

    echo "* nvector_array.ml -> nvector_array.cmo"
    ${OCAMLC} -c nvector_array.ml || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* nvector_array.ml -> nvector_array.cmx"
	${OCAMLOPT} -c ${OCAMLOPTFLAGS} nvector_array.ml || exit 1
    fi

    echo "* ... -> nvector_array.cma"
    ${OCAMLC} -a -o nvector_array.cma -custom nvector_array.cmo || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* ... -> nvector.cmxa"
	${OCAMLOPT} -a -o nvector_array.cmxa ${OCAMLOPTFLAGS} \
	    nvector_array.cmx || exit 1
    fi

    echo "* solvelucy.mli -> solvelucy.cmi"
    ${OCAMLC} solvelucy.mli || exit 1

    echo "* solvelucy.ml -> solvelucy.cmo"
    ${OCAMLC} -c solvelucy.ml || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* solvelucy.ml -> solvelucy.cmx"
	${OCAMLOPT} ${OCAMLOPTFLAGS} -c solvelucy.ml || exit 1
    fi

    # EXAMPLES

    cd examples/

    echo "* examples: showball.mli -> showball.cmi"
    ${OCAMLC} showball.mli || exit 1

    echo "* examples: showball.ml -> showball.cmo"
    ${OCAMLC} -c showball.ml || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* examples: showball.ml -> showball.cmx"
	${OCAMLOPT} -c ${OCAMLOPTFLAGS} showball.ml || exit 1
    fi

    echo "* examples: ball.ml -> ball"
    ${OCAMLC} -o ball -I $LIB -I .. \
	bigarray.cma unix.cma graphics.cma \
	cvode.cma showball.cmo ball.ml || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* examples: ball.ml -> ball.opt"
	${OCAMLOPT} -o ball.opt -I $LIB -I .. ${OCAMLOPTFLAGS} \
	    bigarray.cmxa unix.cmxa graphics.cmxa \
	    cvode.cmxa showball.cmx ball.ml || exit 1
    fi

    for f in $LUCYSOLVE_EXAMPLES; do
	echo "* examples: $f.ml -> $f"
	${OCAMLC} -o $f -I $LIB -I .. \
	    unix.cma bigarray.cma cvode.cma solvelucy.cmo $f.ml || exit 1

	if [ "${OCAMLOPT}" != "" ]; then
	    echo "* examples: $f.ml -> $f.opt"
	    ${OCAMLOPT} -o $f.opt -I $LIB -I .. ${OCAMLOPTFLAGS} \
		unix.cmxa bigarray.cmxa cvode.cmxa solvelucy.cmx $f.ml || exit 1
	fi

	plot_example $f
    done

    for f in $BASIC_EXAMPLES; do
	echo "* examples: $f.ml -> $f"
	${OCAMLC} -o $f -I $LIB -I .. \
	    unix.cma bigarray.cma cvode.cma $f.ml || exit 1

	if [ "${OCAMLOPT}" != "" ]; then
	    echo "* examples: $f.ml -> $f.opt"
	    ${OCAMLOPT} -o $f.opt -I $LIB -I .. ${OCAMLOPTFLAGS} \
		unix.cmxa bigarray.cmxa cvode.cmxa $f.ml || exit 1
	fi

	plot_example $f
    done

    # SUNDIALS EXAMPLES

    cd sundials/

    for f in $SUNDIALS_EXAMPLES; do
	echo "* examples/sundials: $f.ml -> $f"
	${OCAMLC} -o $f -I $LIB -I ../.. \
	    unix.cma bigarray.cma cvode.cma $f.ml \
	    || exit 1

	if [ "${OCAMLOPT}" != "" ]; then
	    echo "* examples/sundials: $f.ml -> $f.opt"
	    ${OCAMLOPT} -o $f.opt -I $LIB -I ../.. ${OCAMLOPTFLAGS} \
		unix.cmxa bigarray.cmxa cvode.cmxa $f.ml \
		|| exit 1
	fi

	plot_example $f
    done

    # SUNDIALS NVECTOR EXAMPLES

    for f in $SUNDIALS_NVECTOR_EXAMPLES; do
	echo "* examples/sundials: $f.ml -> $f"
	${OCAMLC} -o $f -I $LIB -I ../.. \
	    unix.cma bigarray.cma cvode.cma \
	    nvector.cma nvector_array.cma $f.ml \
	    || exit 1

	if [ "${OCAMLOPT}" != "" ]; then
	    echo "* examples/sundials: $f.ml -> $f.opt"
	    ${OCAMLOPT} -o $f.opt -I $LIB -I ../.. ${OCAMLOPTFLAGS} \
		unix.cmxa bigarray.cmxa cvode.cmxa \
		nvector.cmxa nvector_array.cmxa $f.ml \
		|| exit 1
	fi

	plot_example $f
    done

    ;;

esac

