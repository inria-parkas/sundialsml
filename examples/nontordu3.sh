NAME=nontordu3

./${NAME} \
    | sed '/^R/d;s/^[CD]:  //' \
    | head -n 10000 \
    >  ${NAME}.log

cat | gnuplot -persist << END_GNUPLOT
# set terminal postscript color
set title "${NAME}"
set grid
set xtics 1.0
plot "${NAME}.log" using 1:2 with lines title "x",\\
     "${NAME}.log" using 1:3 with lines title "y"
END_GNUPLOT

