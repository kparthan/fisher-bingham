set terminal post eps color enhanced 
box_width=0.20
set style fill solid 0.25 noborder
set style boxplot outliers pointtype 7
set style data boxplot
set boxwidth box_width #relative
set pointsize 0.5
unset key
set border 2 
set xtics nomirror
set ytics nomirror

set xlabel "{/Symbol k}\n"
set ylabel "KL-divergence"
set xlabel font "Times-Roman, 25"
set ylabel font "Times-Roman, 25"
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
#set xtics ("1.10" 1, "1.15" 2, "1.20" 3, "1.25" 4, "1.30" 5, "1.35" 6, "1.40" 7, "1.45" 8, "1.50" 9, "1.55" 10, "1.60" 11) scale 0.0
set xtics ("10" 1, "50" 2, "100" 3) scale 0.0

#set yr[:0.45]
#set label "Proposed" at (1)+0.5,0.44 font "Times-Roman, 25" 
#set arrow from (1)-0.2,0.44 to (1)+0.5,0.44 lw 2 lc rgb "red" nohead 
#set label "FJ" at (1)+0.5,0.41 font "Times-Roman, 25" 
#set arrow from (1)-0.2,0.41 to (1)+0.5,0.41 lw 2 lc rgb "blue" nohead 
d_width=0.5*box_width

set output "N_100_new32_prior/boxplot_kldiv_e_0.1.eps" 

plot  "N_100_new32_prior/k_10_e_0.1/kldivs" using ((1)-3*d_width):1 lt 1 lc rgb "red", \
      "N_100_new32_prior/k_10_e_0.1/kldivs" using ((1)-d_width):2 lt 1 lc rgb "blue", \
      "N_100_new32_prior/k_10_e_0.1/kldivs" using ((1)+d_width):3 lt 1 lc rgb "dark-green", \
      "N_100_new32_prior/k_10_e_0.1/kldivs" using ((1)+3*d_width):4 lt 1 lc rgb "black", \
      "N_100_new32_prior/k_50_e_0.1/kldivs" using ((2)-3*d_width):1 lt 1 lc rgb "red", \
      "N_100_new32_prior/k_50_e_0.1/kldivs" using ((2)-d_width):2 lt 1 lc rgb "blue", \
      "N_100_new32_prior/k_50_e_0.1/kldivs" using ((2)+d_width):3 lt 1 lc rgb "dark-green", \
      "N_100_new32_prior/k_50_e_0.1/kldivs" using ((2)+3*d_width):4 lt 1 lc rgb "black", \
      "N_100_new32_prior/k_100_e_0.1/kldivs" using ((3)-3*d_width):1 lt 1 lc rgb "red", \
      "N_100_new32_prior/k_100_e_0.1/kldivs" using ((3)-d_width):2 lt 1 lc rgb "blue", \
      "N_100_new32_prior/k_100_e_0.1/kldivs" using ((3)+d_width):3 lt 1 lc rgb "dark-green", \
      "N_100_new32_prior/k_100_e_0.1/kldivs" using ((3)+3*d_width):4 lt 1 lc rgb "black"
