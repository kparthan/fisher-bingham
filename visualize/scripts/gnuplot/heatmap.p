#set terminal postscript eps enhanced monochrome 
set terminal postscript eps enhanced color

set output "heatmap.eps" 

set xlabel "x" 
set xrange [0:30] 
set ylabel "y" 
set yrange [0:1]
set zlabel "z" 

unset key

#set ticslevel 0 
#set dgrid3d ,,5
#set hidden3d
set view 42,75 
#set view 27,66 

set pm3d at b 
##set pm3d at b map 
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
splot "../../sampled_data/prior2d_posterior2.dat" using 1:2:3:3 with lines lc rgb "black"

