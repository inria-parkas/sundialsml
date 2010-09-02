
./nontordu2 \
    | sed '/^R/d;s/^[CD]:  //' \
    >  nontordu2.log

cat | gnuplot -persist << END_GNUPLOT
# set terminal postscript color
set title "nontordu2"
set grid
set xtics 1.0
plot "nontordu2.log" using 1:2 with lines title "x",\\
     "nontordu2.log" using 1:3 with lines title "y",\\
     "nontordu2.log" using 1:4 with lines title "z"
END_GNUPLOT

