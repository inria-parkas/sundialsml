#!/bin/sh

envs() {
    echo "You might have to set these environment variables to adjust the output:"
    echo "FONT       - Font name, comma, then size, like: Arial,10"
    echo "SIZE       - Canvas size; see \"help set term size\" in gnuplot."
    echo "BOXCOLOR   - Color of the plot; either a color name or #rrggbb."
    echo "DOTTYPE    - Style (\"pointtype\") of C runtime dots."
    echo "DOTSIZE    - Style (\"pointsize\") of C runtime dots."
    echo "LMARGIN    - How much space (in %) for text at the left."
    echo "RMARGIN    - How much space (in %) for text on the right."
    echo "BMARGIN    - How much space (in %) for text at the bottom."
    echo "TMARGIN    - How much space (in %) for text on the top."
    echo "KEYWIDTH   - Key width increment; see \"help set key\" in gnuplot."
}

if [ "x$1" = x ]; then
    echo "$0 <data file>"
    printf "Makes a plot.  " # in lieu of echo -n
    envs
    echo "Other options:"
    echo "STYLE    - Set this to \"boxplot\" to get a box plot."
    echo "TITLE    - The title."
    echo "TERMINAL - x11, png, jpg, pdf, etc."
    echo "OUTPUT   - Output file, if applicable."
    echo "PAUSE    - Executed at the end; defaults to \"pause mouse\" if TERMINAL=x11."
    echo "YMAX     - Top of the range of left axis."
    echo "Y2MAX    - Top of the range of right axis."
    exit 0
fi

if [ "x$1" = "x--explain-vars" ]; then
    envs
    exit 0
fi

AWK=${AWK:-awk}

YMAX=${YMAX:-*}
Y2MAX=${Y2MAX:-*}

if [ "x$LMARGIN" = x ]; then
    SET_LMARGIN=
else
    LMARGIN=`${AWK} "BEGIN { print ($LMARGIN / 100) }" < /dev/null`
    SET_LMARGIN="set lmargin ${LMARGIN}"
fi

if [ "x$BMARGIN" = x ]; then
    SET_BMARGIN=
else
    BMARGIN=`${AWK} "BEGIN { print ($BMARGIN / 100) }" < /dev/null`
    SET_BMARGIN="set bmargin ${BMARGIN}"
fi

if [ "x$RMARGIN" = x ]; then
    SET_RMARGIN=
else
    RMARGIN=`${AWK} "BEGIN { print (1 - $RMARGIN / 100) }" < /dev/null`
    SET_RMARGIN="set rmargin ${RMARGIN}"
fi

if [ "x$TMARGIN" = x ]; then
    SET_TMARGIN=
else
    TMARGIN=`${AWK} "BEGIN { print (1 - $TMARGIN / 100) }" < /dev/null`
    SET_TMARGIN="set tmargin ${TMARGIN}"
fi

FONT=${FONT:-Arial,7}

case x$FONT in
    x)    SET_FONT=;;
    *)    SET_FONT="set termoption font '$FONT'";;
esac

case x$TERMINAL in
    x)    SET_TERMINAL="set terminal wxt"
          PAUSE=${PAUSE:-"pause mouse"};;
    xwxt) SET_TERMINAL="set terminal wxt"
          PAUSE=${PAUSE:-"pause mouse"};;
    xx11) SET_TERMINAL="set terminal x11"
          PAUSE=${PAUSE:-"pause mouse"};;
    *)    SET_TERMINAL="set terminal $TERMINAL";;
esac

if [ "x$SIZE" != x ]; then
    SET_TERMINAL="$SET_TERMINAL size $SIZE"
fi

KEY_WIDTH=${KEY_WIDTH:--17}
if [ x$KEY_WIDTH != x ]; then
    SET_KEY_WIDTH="set key width $KEY_WIDTH"
fi

if [ "x$OUTPUT" = x ]; then
    SET_OUTPUT=
else
    SET_OUTPUT="set output '$OUTPUT'"
fi

# Compute placement of boxes.  If multiple files are given, cluster
# them.
#
# Given:
# n = number of files
# g = gap between clusters
# c = center of label corresponding to cluster
#
# Compute:
# w = width of each box
# x0 = position of leftmost box's left edge, expressed in cluster's
#      coordinate system (see below)
#
# We have nw+g = 1, so w = (1-g)/n.
#
# In each cluster, the origin of the coordinate system is at c-w/2.
# The left end of the cluster, including half a gap, is at c-1/2.
# Then x0 = c-1/2+g/2.  Expressed in cluster's coordinate system,
# this is c-1/2+g/2 - (c-w/2) = (w+g-1)/2.
GAP=${GAP:-0.5}
w=`${AWK} "BEGIN { printf \"%f\", ((1-${GAP}) / $#) }" < /dev/null`
x0=`${AWK} "BEGIN { printf \"%f\", (($w+${GAP}-1)/2) }" < /dev/null`

if [ "x$POINTSIZE" = x ]; then
    SET_POINTSIZE="set pointsize $w"
else
    SET_POINTSIZE="set pointsize $POINTSIZE"
fi

# Colors of boxes showing time ratios.
if [ "x$BOXCOLORS" = x ]; then
    BOXCOLORS="${BOXCOLOR:-${BOXCOLOR1:-red}} ${BOXCOLOR2:-royalblue} \
               ${BOXCOLOR3:-green} ${BOXCOLOR4:-gold}"
fi

# Colors of dots showing C times.
if [ "x$DOTCOLORS" = x ]; then
    DOTCOLORS="${DOTCOLOR:-${DOTCOLOR1:-black}} ${DOTCOLOR2:-palegreen} \
               ${DOTCOLOR3:-blue} ${DOTCOLOR4:-silver}"
fi

Y2LABEL='C running time / repetition [seconds]'
YLABEL='running time: OCaml / C'

SET_COMMON="$SET_TERMINAL; $SET_OUTPUT; $SET_FONT; set title '$TITLE'"
SET_COMMON="$SET_COMMON; $SET_LMARGIN; $SET_RMARGIN; $SET_BMARGIN"
SET_COMMON="$SET_COMMON; $SET_TMARGIN; $SET_KEY_WIDTH; $SET_POINTSIZE"
SET_COMMON="$SET_COMMON; set ytics nomirror; set y2tics nomirror"
SET_COMMON="$SET_COMMON; set ylabel '${YLABEL}' rotate by -90"
SET_COMMON="$SET_COMMON; set y2label '${Y2LABEL}' rotate by -90"
SET_COMMON="$SET_COMMON; set yrange [0:${YMAX}]; set y2range [0:${Y2MAX}]"
SET_COMMON="$SET_COMMON; set xtics rotate by -90 left"
SET_COMMON="$SET_COMMON; files='$@'"
SET_COMMON="$SET_COMMON; BOXCOLORS='${BOXCOLORS}'; DOTCOLORS='${DOTCOLORS}'"
SET_COMMON="$SET_COMMON; DOTSIZE='${DOTSIZE:-.3}'; DOTTYPE='${DOTTYPE:-7}'"

# C median points must be plotted after boxes, but their key looks
# better above the boxes' key.
SET_COMMON="$SET_COMMON; set key inv"

HORIZ_LINES=
for y in 1 1.5 2; do
    [ "x$HORIZ_LINES" = x ] || HORIZ_LINES="${HORIZ_LINES},"
    HORIZ_LINES="${HORIZ_LINES} $y with lines lc rgb \"#bbbbbb\" notitle"
done

crunch=`echo $0 | sed -e 's#plot\.sh#crunchperf#'`


if [ "x$STYLE" = xboxplot ]; then

LABELCMD="('< $crunch -S $1') u (1):xticlabels(6) lt -3 notitle"

if [ "$#" -eq 1 ]; then
    BOXCMD="'$1' using (0):(\$6):(0.5):7 w boxplot \
               pointtype 2 lc rgb word(BOXCOLORS,1) \
              title 'OCaml time / C time (left axis)'"
    DOTCMD="'$1' using 1:(\$3/\$2) \
               with points pointtype 7 lc rgb 'black' \
               title 'C time / rep (right axis)' axes x1y2"
else
    BOXCMD="for [i=1:words(files)] for [j=0:N-1] \
              word(files,i) index j u (\$1+$x0+(i-1)*$w):(\$6):($w) w boxplot \
              pointtype 2 lc rgb word(BOXCOLORS,i) \
              title word(files,i+j*words(files))"
    DOTCMD="for [i=1:words(files)] \
              ('< $crunch -S '.word(files,i)) \
              u (\$1+($x0)+(i-1)*$w):(\$4/\$2) \
              w points pointtype 7 \
              lc rgb word(DOTCOLORS,i) \
              notitle axes x1y2"
fi

gnuplot <<EOF
$SET_COMMON

# draw vertical lines
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
stats '$1' noout
N=STATS_blocks
set xrange [-$w:STATS_blocks+$w]
set bars $w

# plot the whole set with boxplot, plot each data set's
# median c time / reps with points
plot ${LABELCMD}, ${BOXCMD}, ${DOTCMD}, ${HORIZ_LINES}
${PAUSE}
EOF

else                            # if STYLE != boxplot, do a bar chart

LABELCMD="('< $crunch -S $1') u (1):xticlabels(6) linetype -3 notitle"
if [ "$#" -eq 1 ]; then
    BOXCMD="'< $crunch -S $1' u 1:5 w boxes \
            title 'OCaml time / C time (left axis)' \
            lc rgb word(BOXCOLORS,1)"
    DOTCMD="'< $crunch -S $1' u 1:(\$4/\$2) \
	    w points pointsize DOTSIZE pointtype DOTTYPE \
            lc rgb word(DOTCOLORS,1) \
            title 'C time / rep (right axis)' axes x1y2"
else
    BOXCMD="for [i=1:words(files)] \
              ('< $crunch -S '.word(files,i)) \
              u (\$1+($x0)+(i-1)*$w):5 w boxes \
              title word(files,i) \
              lc rgb word(BOXCOLORS,i)"
    DOTCMD="for [i=1:words(files)] \
              ('< $crunch -S '.word(files,i)) \
              u (\$1+($x0)+(i-1)*$w):(\$4/\$2) \
              w points pointtype 7 \
              lc rgb word(DOTCOLORS,i) \
              notitle axes x1y2"
fi

gnuplot <<EOF
$SET_COMMON

# Bar chart-specific setup
set boxwidth $w
set style fill solid

plot ${LABELCMD}, ${BOXCMD}, ${DOTCMD}, ${HORIZ_LINES}
${PAUSE}
EOF
fi
