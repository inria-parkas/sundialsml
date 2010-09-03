
# todo:
# * use a proper build system
# * also allow for native code
# * also allow for generation and dynamic linking of cvode_serial.so

CC=cc
AR=ar
OCAMLC=ocamlc
LIB=/usr/local/lib
OCAML_INCLUDE=`${OCAMLC} -where`

BASIC_EXAMPLES="discontinuous sincos"
LUCYSOLVE_EXAMPLES="nontordu nontordu2 nontordu3 sincos_lucyf"
SUNDIALS_EXAMPLES="cvRoberts_dns cvAdvDiff_bnd"

case $1 in
clean)
    rm -f cvode_serial.o cvode_serial_bp.o libcvode_serial.a
    rm -f cvode_serial.cmi cvode_serial.cmo
    rm -f solvelucy.cmi solvelucy.cmo
    rm -f cvode_serial.cma

    rm -f examples/ball.cmi examples/ball.cmo
    rm -f examples/showball.cmi examples/showball.cmo examples/showball.cma
    rm -f examples/sincos.cmi examples/sincos.cmo
    rm -f examples/sincos_lucyf.cmi examples/sincos_lucyf.cmo
    rm -f examples/sincos examples/sincos_lucyf examples/ball

    rm -f examples/nontordu examples/nontordu.cmo
    rm -f examples/nontordu.cmi

    for f in $BASIC_EXAMPLES $LUCYSOLVE_EXAMPLES; do
	rm -f examples/$f.cmo;
	rm -f examples/$f.cmi;
	rm -f examples/$f;
    done

    for f in $SUNDIALS_EXAMPLES; do
	rm -f examples/sundials/$f.cmo;
	rm -f examples/sundials/$f.cmi;
	rm -f examples/sundials/$f;
    done
    ;;

*)
    echo "* cvode_serial.c -> cvode_serial.o"
    ${CC} -I $OCAML_INCLUDE -c cvode_serial.c || exit 1

    echo "* cvode_serial_bp.c -> cvode_serial_bp.o"
    ${CC} -I $OCAML_INCLUDE -c cvode_serial_bp.c || exit 1

    echo "* cvode_serial.o -> libcvode_serial.a"
    ${AR} rc libcvode_serial.a cvode_serial.o cvode_serial_bp.o || exit 1

    echo "* cvode_serial.mli -> cvode_serial.cmi"
    ${OCAMLC} cvode_serial.mli || exit 1

    echo "* cvode_serial.ml -> cvode_serial.cmo"
    ${OCAMLC} -c cvode_serial.ml || exit 1

    echo "* ... -> cvode_serial.cma"
    ${OCAMLC} -a -o cvode_serial.cma -custom cvode_serial.cmo \
	-cclib -lsundials_cvode \
	-cclib -lsundials_nvecserial \
	-cclib -lcvode_serial || exit 1

    echo "* solvelucy.mli -> solvelucy.cmi"
    ${OCAMLC} solvelucy.mli || exit 1

    echo "* solvelucy.ml -> solvelucy.cmo"
    ${OCAMLC} -c solvelucy.ml || exit 1

    # EXAMPLES

    cd examples/

    echo "* examples: showball.mli -> showball.cmi"
    ${OCAMLC} showball.mli || exit 1

    echo "* examples: showball.ml -> showball.cmo"
    ${OCAMLC} -c showball.ml || exit 1

    echo "* examples: ... -> showball.cma"
    ${OCAMLC} -a -o showball.cma unix.cma graphics.cma showball.cmo || exit 1

    echo "* examples: ball.ml -> ball"
    ${OCAMLC} -o ball -I $LIB -I .. \
	bigarray.cma unix.cma \
	cvode_serial.cma showball.cma ball.ml || exit 1

    for f in $LUCYSOLVE_EXAMPLES; do
	echo "* examples: $f.ml -> $f"
	${OCAMLC} -o $f -I $LIB -I .. \
	    unix.cma bigarray.cma cvode_serial.cma solvelucy.cmo $f.ml || exit 1
    done

    for f in $BASIC_EXAMPLES; do
	echo "* examples: $f.ml -> $f"
	${OCAMLC} -o $f -I $LIB -I .. \
	    unix.cma bigarray.cma cvode_serial.cma $f.ml || exit 1
    done

    # SUNDIALS EXAMPLES

    cd sundials/

    for f in $SUNDIALS_EXAMPLES; do
	echo "* examples/sundials: $f.ml -> $f"
	${OCAMLC} -o $f -I $LIB -I ../.. \
	    unix.cma bigarray.cma cvode_serial.cma $f.ml || exit 1
    done

    ;;

esac

