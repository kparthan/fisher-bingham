set terminal postscript eps enhanced color

set output "./N_10_prior2/fixed_kappa/kappa_1/errors/biassq_beta.eps"

set style data linespoints
set style fill solid 1.0 noborder
set xlabel "eccentricity"
set ylabel "Bias squared"
plot "./N_10_prior2/fixed_kappa/kappa_1/errors/biassq_beta" using 1:2 t "MOM" lc rgb "red", \
"" using 1:7 t "MAP3 = MLE" lc rgb "blue", \
"" using 1:4 t "MAP1" lc rgb "dark-green", \
"" using 1:5 t "MML" lc rgb "dark-magenta", \
"" using 1:6 t "MAP2" lc rgb "black"
