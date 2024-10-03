set xlabel 'm'
set ylabel 'Error (Pa)'
set output 'P3-23-24-fig1.png'
set term png
set xzeroaxis
set yzeroaxis
set key left top
set key box vertical width 1 height 1 maxcols 1 spacing 1
# set yrange [-0.5:1]
set logscale y

set grid

plot "builds/resE1.dat" index 0 using 1:3 title "Trapezis" with linespoints, \
     "builds/resE1.dat" index 0 using 1:5 title "Euler-McLaurin" with linespoints
