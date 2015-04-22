#set terminal postscript eps enhanced monochrome 
set terminal postscript eps color enhanced 
#set terminal pngcairo size 800,600 enhanced 

set output "heatmap.eps" 

set xlabel "{/Symbol k}" 
set ylabel "{/Symbol b}" 
#set xlabel "z_4" 
#set ylabel "z_5" 
#set zlabel "z" 

set xrange [0:30] 
set yrange [0:15]
set ytics 3

set lmargin 0
set rmargin 0
set tmargin 0 
set bmargin 0

set xlabel font "Times-Roman, 25"
set ylabel font "Times-Roman, 25"
set xtics font "Times-Roman, 20"
set ytics font "Times-Roman, 20"
set xtics nomirror
set ytics nomirror

#set format z "%.1f"
#set label "x 10^{-4}" at 0,-0.01,5

unset colorbox
set border 895

#set ticslevel 0 
#set dgrid3d ,,5
#set hidden3d
set view 42,75 
#set view 27,66 

set pm3d at b 
#set pm3d at b map 
#set contour
##set contour base
#set nosurface
#set isosample 40,40
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")
#splot "../../sampled_data/prior2d_posterior2.dat" using 1:2:($3*1e4) notitle with lines lc rgb "black"
splot "../../sampled_data/prior2d_posterior1.dat" using 1:2:3 notitle with lines lc rgb "black"
