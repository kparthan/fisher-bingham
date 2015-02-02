set terminal post eps color enhanced
set autoscale	# scale axes automatically
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set style data linespoints

set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set xlabel font "Times-Roman, 25"
set ylabel font "Times-Roman, 25"

set xlabel "eccentricity\n"
set ylabel "MSE (kappa)"

set xr[0.1:]
set yr[:]

set xtics 0.1

#set key at 1.55,30
#set key spacing 3.0 font "Times-Roman, 25"

set output "mse_k_10.eps"

plot "test.dat" using 1:2 title "MOMENT" lc rgb "red", \
     "" using 1:3 title "MLE" lc rgb "blue", \
     "" using 1:4 title "MAP" lc rgb "dark-green", \
     "" using 1:5 title "MML" lc rgb "black"
