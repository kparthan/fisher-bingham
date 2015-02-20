set terminal post eps color enhanced 

box_width=0.10
set style fill solid 0.25 noborder
set grid y
set style data histograms
set boxwidth box_width 

set xlabel "eccentricity\n"
set ylabel "KL-divergence"
set xlabel font "Times-Roman, 25"
set ylabel font "Times-Roman, 25"
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"

#set xtics ("0.1" 1, "0.5" 2, "0.9" 3) scale 0.0
#set xtics ("0.1" 1, "0.5" 2) scale 0.0
#set xtics ("0.1" 1) scale 0.0
set xtics 0.1,0.4
d_width=0.5*box_width

set arrow from -0.2,0 to 1.2,0 nohead

set output "N_25_new32_prior/boxplot_kldivs_diff_k_10.eps" 

plot "N_25_new32_prior/kldivs_diff_k_10.dat" using ($1-2*d_width):2 with boxes lc rgb "red" notitle, \
     "" using ($1):3 with boxes lc rgb "blue" notitle, \
     "" using ($1+2*d_width):4 with boxes lc rgb "dark-green" notitle, \
     "" using ($1-2*d_width):5 with boxes lc rgb "red" notitle, \
     "" using ($1):6 with boxes lc rgb "blue" notitle, \
     "" using ($1+2*d_width):7 with boxes lc rgb "dark-green" notitle
