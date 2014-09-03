#!/bin/sh

envs() {
    echo "You might have to set these environment variables to adjust the output:"
    echo "FONT       - Font name, comma, then size, like: Arial,10"
    echo "TITLE_FONT - Like FONT, but only used for the title."
    echo "SIZE       - Canvas size; see \"help set term size\" in gnuplot."
    echo "LMARGIN    - How much space (in %) for text at the left."
    echo "RMARGIN    - How much space (in %) for text on the right."
    echo "BMARGIN    - How much space (in %) for text at the bottom."
    echo "TMARGIN    - How much space (in %) for text on the top."
}

if test "x$1" = x; then
    echo "$0 <data file>"
    printf "Makes a plot.  " # in lieu of echo -n
    envs
    echo "Other options:"
    echo "TITLE    - The title."
    echo "TERMINAL - x11, png, jpg, pdf, etc."
    echo "OUTPUT   - Output file, if applicable."
    echo "PAUSE    - Executed at the end; defaults to \"pause mouse\" if TERMINAL=x11."
    echo "YMAX     - Top of the range of left axis."
    exit 0
fi

if test "x$1" = "x--explain-vars"; then
    envs
    exit 0
fi

AWK=${AWK:-awk}

YMAX=${YMAX:-*}
Y2MAX=${Y2MAX:-*}

if test "x$LMARGIN" = x; then
    LMARGIN=0.1
else
    LMARGIN=`${AWK} "BEGIN { print ($LMARGIN / 100) }"`
fi

if test "x$BMARGIN" = x; then
    BMARGIN=0.33
else
    BMARGIN=`${AWK} "BEGIN { print ($BMARGIN / 100) }"`
fi

if test "x$RMARGIN" = x; then
    RMARGIN=0.93
else
    RMARGIN=`${AWK} "BEGIN { print (1 - $RMARGIN / 100) }"`
fi

if test "x$TMARGIN" = x; then
    TMARGIN=0.93
else
    TMARGIN=`${AWK} "BEGIN { print (1 - $TMARGIN / 100) }"`
fi

TITLE_FONT=${TITLE_FONT:-Arial,8}
FONT=${FONT:-Arial,6}

case x$TITLE_FONT in
    x)    SET_TITLE_FONT=;;
    *)    SET_TITLE_FONT="set title font '$TITLE_FONT'"
esac

case x$FONT in
    x)    SET_FONT=;;
    *)    SET_FONT="set xtics font '$FONT'; \
                    set ytics font '$FONT'; \
                    set y2tics font '$FONT'; \
                    set key font '$FONT'"
esac

case x$TERMINAL in
    x)    SET_TERMINAL="set terminal wxt"
          PAUSE="pause mouse";;
    xwxt) SET_TERMINAL="set terminal wxt"
          PAUSE="pause mouse";;
    xx11) SET_TERMINAL="set terminal x11"
          PAUSE="pause mouse";;
    *)    SET_TERMINAL="set terminal $TERMINAL"
esac

if test "x$SIZE" != x; then
    SET_TERMINAL="$SET_TERMINAL size $SIZE"
fi

KEY_WIDTH=${KEY_WIDTH:--17}
if test x$KEY_WIDTH != x; then
    SET_KEY_WIDTH="set key width $KEY_WIDTH"
fi

if test "x$OUTPUT" = x; then
    SET_OUTPUT=
else
    SET_OUTPUT="set output '$OUTPUT'"
fi

if test "x$POINTSIZE" = x; then
    SET_POINTSIZE=
else
    SET_POINTSIZE="set pointsize $POINTSIZE"
fi

if test "x$STYLE" = xboxplot; then
gnuplot <<EOF
$SET_TERMINAL
$SET_OUTPUT
$SET_FONT
$SET_TITLE_FONT
set title '$TITLE'
set lmargin screen ${LMARGIN};
set rmargin screen ${RMARGIN};
set bmargin screen ${BMARGIN};
set tmargin screen ${TMARGIN};
set key inv
$SET_KEY_WIDTH
$SET_POINTSIZE
set xtics rotate by 60 right;
set ytics nomirror;
set y2label 'seconds' font '${FONT}'
set y2tics nomirror;
set yrange [0:${YMAX}];
set y2range [0:${Y2MAX}];

# N = number of data sets
stats "perf.opt.log" noout
N=STATS_blocks

# Draw vertical lines
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"

# plot the whole set with boxplot, plot each data set's
# median C time / reps with points
plot "perf.opt.log" using (0):(\$3/\$4):(0.5):6 with boxplot pointtype 2 lw 0.5 \
       title 'OCaml time / C time (left axis)', \
     "perf.opt.log" index 0 using (0):(\$2/\$1) \
       with points pointtype 3 lw 0.5 lc rgb 'black' \
       title 'C time / rep (right axis)' axes x1y2, \
     for [IDX=1:N-1] "perf.opt.log" index IDX using (IDX):(\$2/\$1) \
       with points pointtype 3 lw 0.5 lc rgb 'black' notitle axes x1y2, \
     1 with lines lc rgb '#bbbbbb' notitle, \
     2 with lines lc rgb '#bbbbbb' notitle
${PAUSE}
EOF
else                            # if STYLE != boxplot, do a bar chart
gnuplot <<EOF
$SET_TERMINAL
$SET_OUTPUT
$SET_FONT
$SET_TITLE_FONT
set title '$TITLE'
set lmargin screen ${LMARGIN};
set rmargin screen ${RMARGIN};
set bmargin screen ${BMARGIN};
set tmargin screen ${TMARGIN};
set key inv
# $SET_KEY_WIDTH
# $SET_POINTSIZE
set xtics rotate by 60 right;
set ytics nomirror;
set y2label 'seconds' font '${FONT}'
set y2tics nomirror;
set yrange [0:${YMAX}];
set y2range [0:${Y2MAX}];

# Bar chart-specific setup
set boxwidth 0.5
set style fill solid

# plot the whole set with boxplot, plot each data set's
# median C time / reps with points
# The grep is needed to merge multiple data sets into one.
plot "< grep '.' $1" using 5:xtic(6) with boxes lc rgb 'red' \
       title 'OCaml time / C time (left axis)', \
     "< grep '.' $1" using (\$2/\$1) with points pointtype 3 lc rgb 'black' \
       title 'C time / rep (right axis)' axes x1y2, \
     1 with lines lc rgb "#bbbbbb" notitle, \
     2 with lines lc rgb "#bbbbbb" notitle
$PAUSE
EOF
fi
