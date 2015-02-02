set terminal post eps color enhanced 
box_width=0.25
set style fill solid 0.25 noborder
set style boxplot outliers pointtype 7
set style data boxplot
set boxwidth box_width #relative
set pointsize 0.5
unset key
set border 2 
set xtics nomirror
set ytics nomirror

set xlabel "Separation {/Symbol d}\n"
set ylabel "Empirical KL-divergence"
set xlabel font "Times-Roman, 25"
set ylabel font "Times-Roman, 25"
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set xtics ("1.10" 1, "1.15" 2, "1.20" 3, "1.25" 4, "1.30" 5, "1.35" 6, "1.40" 7, "1.45" 8, "1.50" 9, "1.55" 10, "1.60" 11) scale 0.0
#set xtics ("1.10" 1, "1.15" 2, "1.20" 3, "1.25" 4, "1.30" 5, "1.35" 6, "1.40" 7) scale 0.0

#set yr[:0.45]
#set label "Proposed" at (1)+0.5,0.44 font "Times-Roman, 25" 
#set arrow from (1)-0.2,0.44 to (1)+0.5,0.44 lw 2 lc rgb "red" nohead 
#set label "FJ" at (1)+0.5,0.41 font "Times-Roman, 25" 
#set arrow from (1)-0.2,0.41 to (1)+0.5,0.41 lw 2 lc rgb "blue" nohead 
d_width=0.5*box_width

set output "comparisons_boxplot_kldiv_10d_exp2a.eps" 
plot  "kldivs2_delta_1.10" using ((1)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.10" using ((1)+d_width):5 lt 1 lc rgb "blue", \
      "kldivs2_delta_1.15" using ((2)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.15" using ((2)+d_width):5 lt 1 lc rgb "blue", \
      "kldivs2_delta_1.20" using ((3)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.20" using ((3)+d_width):5 lt 1 lc rgb "blue", \
      "kldivs2_delta_1.25" using ((4)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.25" using ((4)+d_width):5 lt 1 lc rgb "blue", \
      "kldivs2_delta_1.30" using ((5)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.30" using ((5)+d_width):5 lt 1 lc rgb "blue", \
      "kldivs2_delta_1.35" using ((6)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.35" using ((6)+d_width):5 lt 1 lc rgb "blue", \
      "kldivs2_delta_1.40" using ((7)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.40" using ((7)+d_width):5 lt 1 lc rgb "blue", \
      "kldivs2_delta_1.45" using ((8)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.45" using ((8)+d_width):5 lt 1 lc rgb "blue", \
      "kldivs2_delta_1.50" using ((9)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.50" using ((9)+d_width):5 lt 1 lc rgb "blue", \
      "kldivs2_delta_1.55" using ((10)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.55" using ((10)+d_width):5 lt 1 lc rgb "blue", \
      "kldivs2_delta_1.60" using ((11)-d_width):1 lt 1 lc rgb "red", \
      "kldivs2_delta_1.60" using ((11)+d_width):5 lt 1 lc rgb "blue"
