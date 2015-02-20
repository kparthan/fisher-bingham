set terminal post eps color enhanced 

set xlabel "eccentricity\n"
set xlabel font "Times-Roman, 25"
set ylabel font "Times-Roman, 25"
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"

set style data linespoints

set key top left
set ylabel "Avg. KL-divergence"
set output "N_1000_new32_prior/k_10_avg_kldivs.eps" 
plot "N_1000_new32_prior/k_10_avg_kldivs.dat" using  2:3 title "Moment" lc rgb "red", \
     "" using  2:4 title "MLE" lc rgb "blue", \
     "" using  2:5 title "MAP" lc rgb "dark-green", \
     "" using  2:6 title "MML" lc rgb "black"

set key top right
set ylabel "% KL-divergence Wins"
set output "N_1000_new32_prior/k_10_wins_kldivs.eps" 
plot "N_1000_new32_prior/k_10_wins_kldivs.dat" using  2:3 title "Moment" lc rgb "red", \
     "" using  2:4 title "MLE" lc rgb "blue", \
     "" using  2:5 title "MAP" lc rgb "dark-green", \
     "" using  2:6 title "MML" lc rgb "black"
