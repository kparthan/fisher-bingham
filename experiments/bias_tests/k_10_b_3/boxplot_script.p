set terminal postscript eps enhanced
#set style fill solid 0.25 border -1
set style boxplot outliers pointtype 7
set style data boxplot
set boxwidth  0.05
set pointsize 0.5
unset key
set border 2

#set xtics nomirror
set ytics nomirror

set xtics ("50" 1, "100" 2, "250" 3, "500" 4, "1000" 5) scale 0.0
set xtics font "Times-Roman, 25"
set ytics font "Times-Roman, 20"

set xlabel "Sample size"
set ylabel "Kappa estimates"
set ylabel font "Times-Roman, 25"

#set title 

set xrange [0:*]
#set yrange [:300]

# plotting kappa

set output "moment_kappa.eps" 
plot "n_50_k_10" using (1):2 lt 1 lc rgb "blue" pointtype 7, \
     "n_100_k_10" using (2):2 lt 1 lc rgb "blue" pointtype 7, \
     "n_250_k_10" using (3):2 lt 1 lc rgb "blue" pointtype 7, \
     "n_500_k_10" using (4):2 lt 1 lc rgb "blue" pointtype 7, \
     "n_1000_k_10" using (5):2 lt 1 lc rgb "blue" pointtype 7, \
     10 with lines lt 4, \
     "median_kappa"  using ($0+1):2 with linespoints lt 4 lc rgb "red"

set output "mle_kappa.eps" 
plot "n_50_k_10" using (1):3 lt 1 lc rgb "blue" pointtype 7, \
     "n_100_k_10" using (2):3 lt 1 lc rgb "blue" pointtype 7, \
     "n_250_k_10" using (3):3 lt 1 lc rgb "blue" pointtype 7, \
     "n_500_k_10" using (4):3 lt 1 lc rgb "blue" pointtype 7, \
     "n_1000_k_10" using (5):3 lt 1 lc rgb "blue" pointtype 7, \
     10 with lines lt 4, \
     "median_kappa"  using ($0+1):3 with linespoints lt 4 lc rgb "red"

set output "map_kappa.eps" 
plot "n_50_k_10" using (1):4 lt 1 lc rgb "blue" pointtype 7, \
     "n_100_k_10" using (2):4 lt 1 lc rgb "blue" pointtype 7, \
     "n_250_k_10" using (3):4 lt 1 lc rgb "blue" pointtype 7, \
     "n_500_k_10" using (4):4 lt 1 lc rgb "blue" pointtype 7, \
     "n_1000_k_10" using (5):4 lt 1 lc rgb "blue" pointtype 7, \
     10 with lines lt 4, \
     "median_kappa"  using ($0+1):4 with linespoints lt 4 lc rgb "red"

set output "mml_kappa.eps" 
plot "n_50_k_10" using (1):5 lt 1 lc rgb "blue" pointtype 7, \
     "n_100_k_10" using (2):5 lt 1 lc rgb "blue" pointtype 7, \
     "n_250_k_10" using (3):5 lt 1 lc rgb "blue" pointtype 7, \
     "n_500_k_10" using (4):5 lt 1 lc rgb "blue" pointtype 7, \
     "n_1000_k_10" using (5):5 lt 1 lc rgb "blue" pointtype 7, \
     10 with lines lt 4, \
     "median_kappa"  using ($0+1):5 with linespoints lt 4 lc rgb "red"

# plotting beta

set output "moment_beta.eps" 
plot "n_50_b_3" using (1):2 lt 1 lc rgb "blue" pointtype 7, \
     "n_100_b_3" using (2):2 lt 1 lc rgb "blue" pointtype 7, \
     "n_250_b_3" using (3):2 lt 1 lc rgb "blue" pointtype 7, \
     "n_500_b_3" using (4):2 lt 1 lc rgb "blue" pointtype 7, \
     "n_1000_b_3" using (5):2 lt 1 lc rgb "blue" pointtype 7, \
     3 with lines lt 4, \
     "median_beta"  using ($0+1):2 with linespoints lt 4 lc rgb "red"

set output "mle_beta.eps" 
plot "n_50_b_3" using (1):3 lt 1 lc rgb "blue" pointtype 7, \
     "n_100_b_3" using (2):3 lt 1 lc rgb "blue" pointtype 7, \
     "n_250_b_3" using (3):3 lt 1 lc rgb "blue" pointtype 7, \
     "n_500_b_3" using (4):3 lt 1 lc rgb "blue" pointtype 7, \
     "n_1000_b_3" using (5):3 lt 1 lc rgb "blue" pointtype 7, \
     3 with lines lt 4, \
     "median_beta"  using ($0+1):3 with linespoints lt 4 lc rgb "red"

set output "map_beta.eps" 
plot "n_50_b_3" using (1):4 lt 1 lc rgb "blue" pointtype 7, \
     "n_100_b_3" using (2):4 lt 1 lc rgb "blue" pointtype 7, \
     "n_250_b_3" using (3):4 lt 1 lc rgb "blue" pointtype 7, \
     "n_500_b_3" using (4):4 lt 1 lc rgb "blue" pointtype 7, \
     "n_1000_b_3" using (5):4 lt 1 lc rgb "blue" pointtype 7, \
     3 with lines lt 4, \
     "median_beta"  using ($0+1):4 with linespoints lt 4 lc rgb "red"

set output "mml_beta.eps" 
plot "n_50_b_3" using (1):5 lt 1 lc rgb "blue" pointtype 7, \
     "n_100_b_3" using (2):5 lt 1 lc rgb "blue" pointtype 7, \
     "n_250_b_3" using (3):5 lt 1 lc rgb "blue" pointtype 7, \
     "n_500_b_3" using (4):5 lt 1 lc rgb "blue" pointtype 7, \
     "n_1000_b_3" using (5):5 lt 1 lc rgb "blue" pointtype 7, \
     3 with lines lt 4, \
     "median_beta"  using ($0+1):5 with linespoints lt 4 lc rgb "red"
