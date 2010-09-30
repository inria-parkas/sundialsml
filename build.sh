
# todo:
# * use a proper build system
# * also allow for native code
# * also allow for generation and dynamic linking of cvode_serial.so

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
LUCYSOLVE_EXAMPLES="nontordu nontordu2 nontordu3 sincos_lucyf billiard1d"
SUNDIALS_EXAMPLES="cvRoberts_dns cvAdvDiff_bnd"
PLOT="examples/billiard1d examples/nontordu2 examples/nontordu3 "

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
    rm -f cvode_serial.o cvode_serial_bp.o libcvode_serial.a
    rm -f cvode_serial.cmi cvode_serial.cmo
    rm -f solvelucy.cmi solvelucy.cmo
    rm -f cvode_serial.cma

    rm -f cvode_serial.cmx cvode_serial.cmxa cvode_serial.a
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
    echo "* cvode_serial.c -> cvode_serial.o"
    ${CC} -I $OCAML_INCLUDE ${CFLAGS} -c cvode_serial.c || exit 1

    echo "* cvode_serial_bp.c -> cvode_serial_bp.o"
    ${CC} -I $OCAML_INCLUDE ${CFLAGS} -c cvode_serial_bp.c || exit 1

    echo "* cvode_serial.o -> libcvode_serial.a"
    ${AR} rc libcvode_serial.a cvode_serial.o cvode_serial_bp.o || exit 1

    echo "* cvode_serial.mli -> cvode_serial.cmi"
    ${OCAMLC} cvode_serial.mli || exit 1

    echo "* cvode_serial.ml -> cvode_serial.cmo"
    ${OCAMLC} -c cvode_serial.ml || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* cvode_serial.ml -> cvode_serial.cmx"
	${OCAMLOPT} -c ${OCAMLOPTFLAGS} cvode_serial.ml || exit 1
    fi

    echo "* ... -> cvode_serial.cma"
    ${OCAMLC} -a -o cvode_serial.cma -custom cvode_serial.cmo \
	-cclib -lsundials_cvode \
	-cclib -lsundials_nvecserial \
	${LAPACK_LIB} \
	-cclib -lcvode_serial || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* ... -> cvode_serial.cmxa"
	${OCAMLOPT} -a -o cvode_serial.cmxa ${OCAMLOPTFLAGS} \
	    cvode_serial.cmx \
	    -cclib -lsundials_cvode \
	    -cclib -lsundials_nvecserial \
	    ${LAPACK_LIB} \
	    -cclib -lcvode_serial || exit 1
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
	cvode_serial.cma showball.cmo ball.ml || exit 1

    if [ "${OCAMLOPT}" != "" ]; then
	echo "* examples: ball.ml -> ball.opt"
	${OCAMLOPT} -o ball.opt -I $LIB -I .. ${OCAMLOPTFLAGS} \
	    bigarray.cmxa unix.cmxa graphics.cmxa \
	    cvode_serial.cmxa showball.cmx ball.ml || exit 1
    fi

    for f in $LUCYSOLVE_EXAMPLES; do
	echo "* examples: $f.ml -> $f"
	${OCAMLC} -o $f -I $LIB -I .. \
	    unix.cma bigarray.cma cvode_serial.cma solvelucy.cmo $f.ml || exit 1

	if [ "${OCAMLOPT}" != "" ]; then
	    echo "* examples: $f.ml -> $f.opt"
	    ${OCAMLOPT} -o $f.opt -I $LIB -I .. ${OCAMLOPTFLAGS} \
		unix.cmxa bigarray.cmxa cvode_serial.cmxa solvelucy.cmx $f.ml || exit 1
	fi

	plot_example $f
    done

    for f in $BASIC_EXAMPLES; do
	echo "* examples: $f.ml -> $f"
	${OCAMLC} -o $f -I $LIB -I .. \
	    unix.cma bigarray.cma cvode_serial.cma $f.ml || exit 1

	if [ "${OCAMLOPT}" != "" ]; then
	    echo "* examples: $f.ml -> $f.opt"
	    ${OCAMLOPT} -o $f.opt -I $LIB -I .. ${OCAMLOPTFLAGS} \
		unix.cmxa bigarray.cmxa cvode_serial.cmxa $f.ml || exit 1
	fi

	plot_example $f
    done

    # SUNDIALS EXAMPLES

    cd sundials/

    for f in $SUNDIALS_EXAMPLES; do
	echo "* examples/sundials: $f.ml -> $f"
	${OCAMLC} -o $f -I $LIB -I ../.. \
	    unix.cma bigarray.cma cvode_serial.cma $f.ml || exit 1

	if [ "${OCAMLOPT}" != "" ]; then
	    echo "* examples/sundials: $f.ml -> $f.opt"
	    ${OCAMLOPT} -o $f.opt -I $LIB -I ../.. ${OCAMLOPTFLAGS} \
		unix.cmxa bigarray.cmxa cvode_serial.cmxa $f.ml || exit 1
	fi

	plot_example $f
    done

    ;;

esac

