
# todo:
# * use a proper build system
# * also allow for native code
# * also allow for generation and dynamic linking of cvode_serial.so

case $1 in
clean)
    rm -f cvode_serial.o cvode_serial_bp.o libcvode_serial.a
    rm -f cvode_serial.cmi cvode_serial.cmo
    rm -f ball.cmi ball.cmo
    rm -f showball.cmi showball.cmo showball.cma
    rm -f solvelucy.cmi solvelucy.cmo
    rm -f cvode_serial.cma
    rm -f sincos.cmi sincos.cmo
    rm -f sincos_lucyf.cmi sincos_lucyf.cmo
    rm -f sincos sincos_lucyf ball
    ;;

*)
    echo "* cvode_serial.c -> cvode_serial.o"
    cc -c cvode_serial.c || exit 1

    echo "* cvode_serial_bp.c -> cvode_serial_bp.o"
    cc -c cvode_serial_bp.c || exit 1

    echo "* cvode_serial.o -> libcvode_serial.a"
    ar rc libcvode_serial.a cvode_serial.o cvode_serial_bp.o || exit 1

    echo "* cvode_serial.mli -> cvode_serial.cmi"
    ocamlc cvode_serial.mli || exit 1

    echo "* cvode_serial.ml -> cvode_serial.cmo"
    ocamlc -c cvode_serial.ml || exit 1

    echo "* ... -> cvode_serial.cma"
    ocamlc -a -o cvode_serial.cma -custom cvode_serial.cmo \
	-cclib -lsundials_cvode \
	-cclib -lsundials_nvecserial \
	-cclib -lcvode_serial || exit 1

    echo "* solvelucy.mli -> solvelucy.cmi"
    ocamlc solvelucy.mli || exit 1

    echo "* solvelucy.ml -> solvelucy.cmo"
    ocamlc -c solvelucy.ml || exit 1

    echo "* sincos.ml -> sincos"
    ocamlc -o sincos -I /usr/local/lib -I . \
	unix.cma bigarray.cma cvode_serial.cma sincos.ml || exit 1

    echo "* sincos_lucyf.ml -> sincos_lucyf"
    ocamlc -o sincos_lucyf -I /usr/local/lib -I . \
	unix.cma bigarray.cma cvode_serial.cma solvelucy.cmo sincos_lucyf.ml || exit 1

    echo "* showball.mli -> showball.cmi"
    ocamlc showball.mli || exit 1

    echo "* showball.ml -> showball.cmo"
    ocamlc -c showball.ml || exit 1

    echo "* ... -> showball.cma"
    ocamlc -a -o showball.cma unix.cma graphics.cma showball.cmo || exit 1

    echo "* ball.ml -> ball"
    ocamlc -o ball -I /usr/local/lib -I . \
	bigarray.cma unix.cma \
	cvode_serial.cma showball.cma ball.ml || exit 1
    ;;

esac

