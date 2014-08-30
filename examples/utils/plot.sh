#!/bin/sh

envs() {
    echo "You might have to set these environment variables to adjust the output:"
    echo "FONT       - Font name, comma, then size, like: Arial,10"
    echo "TITLE_FONT - Like FONT, but only used for the title."
    echo "SIZE       - Canvas size; see \"help size\" in gnuplot."
    echo "BMARGIN    - How much space (in %) for text at the bottom."
    echo "RMARGIN    - How much space (in %) for text on the right."
}

if test "$#" -eq 0; then
    echo "$0 <data file>"
    printf "Makes a plot.  " # in lieu of echo -n
    envs
    echo "Other options:"
    echo "TITLE    - The title."
    echo "TERMINAL - x11, png, jpg, pdf, etc."
    echo "OUTPUT   - Output file, if applicable."
    echo "PAUSE    - Executed at the end; defaults to \"pause mouse\" if TERMINAL=x11."
    exit 0
fi

if test "x$1" = "x--explain-vars"; then
    envs
    exit 0
fi

AWK=${AWK:-awk}

if test "x$RMARGIN" = x; then
    RMARGIN=0.93
else
    RMARGIN=`$(AWK) 'BEGIN { print (1 - $RMARGIN / 100) }'`
fi

if test "x$BMARGIN" = x; then
    BMARGIN=0.33
else
    BMARGIN=`$(AWK) 'BEGIN { print ($BMARGIN / 100) }'`
fi

TITLE_FONT=${TITLE_FONT:-Arial,8}
FONT=${FONT:-Arial,6}

case x$TERMINAL in
    x)    SET_TERMINAL=
          PAUSE="pause mouse";;
    xx11) SET_TERMINAL=
          PAUSE="pause mouse";;
    *)    SET_TERMINAL="set terminal $TERMINAL"
esac

if test "x$OUTPUT" = x; then
    SET_OUTPUT=
else
    SET_OUTPUT="set output '$OUTPUT'"
fi

if test "x$SIZE" = x; then
    SET_SIZE=
else
    SET_SIZE="set size $SIZE"
fi

gnuplot <<EOF
$SET_TERMINAL
$SET_OUTPUT
$SET_SIZE
set title '$TITLE'
set title font '${TITLE_FONT}';
set xtics font '${FONT}';
set ytics font '${FONT}';
set key font '${FONT}';
set key ins vert left top;
set key invert;
set y2tics font '${FONT}';
set bmargin screen ${BMARGIN};
set rmargin screen ${RMARGIN};
set boxwidth 0.5;
set style fill solid;
set yrange [0:*];
set xtics rotate by 60 right;
set logscale y2;
set ytics nomirror;
set y2tics nomirror;
plot "$1" using 4:xtic(5) with boxes lc rgb 'red' \
    title 'OCaml time / C time (left axis)', \
    "$1" using 1 with points \
    pointtype 3 lc rgb 'black' \
    title 'number of repetitions (right axis, log scale)' \
    axes x1y2; $PAUSE
EOF
