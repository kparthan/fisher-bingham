set terminal postscript eps enhanced color
 
set output "stacked_bars.eps"

set key invert reverse left top
 
set grid y
set style data histograms
set style histogram rowstacked
set boxwidth 0.5
set style fill solid 1.0 border -1
set ytics 10 nomirror
set yrange [:100]
set ylabel "% of wins "
set ytics 10
 
plot "test_values.dat" using 2 t "MOM", \
     "" using 3 t "MLE", \
     "" using 4 t "MAP", \
     "" using 5:xtic(1) t "MML"

