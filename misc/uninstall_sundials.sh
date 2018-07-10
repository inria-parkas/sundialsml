#!/bin/sh

PREFIX="${PREFIX:-/usr/local}"
SOLVERS="nvector sunmatrix sunlinsol arkode cvode cvodes ida idas kinsol"

DO_RM="rm -f"
RM="echo ${DO_RM}"
DRYRUN=1
STOP=0

while [ $# -gt 0 ]; do
    case "$1" in
    -h|--help)
	printf "Remove the files installed by Sundials 2.x/3.x.\n"
	printf "By default, simply show the files to be deleted.\n"
	printf "Run with -y to really delete the files.\n"
	printf "Set the install prefix with PREFIX or --prefix <path>\n"
	STOP=1
    ;;
    -y)
	RM="${DO_RM}"
	DRYRUN=0
    ;;
    --prefix)
	PREFIX="$2"
	shift
    ;;
    esac
    shift
done

if [ "$STOP" -eq 1 ]; then
    exit 0
fi

$RM -r "${PREFIX}/include/sundials"

for s in $SOLVERS
do
    $RM -r "${PREFIX}/include/$s"
    $RM -r "${PREFIX}/examples/$s"	# standard examples directory
done

for l in ${PREFIX}/lib/libsundials_*
do
    $RM "$l"
done

if [ "$DRYRUN" -eq 1 ]; then
    printf "** DRY RUN ONLY: rerun with -y to really delete these files.\n" >&2
fi

