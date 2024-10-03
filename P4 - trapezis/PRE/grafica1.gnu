set ylabel 'Divergencia area'
set xlabel 'h'
# set xzeroaxis
set key left top
set key box vertical width 1 height 1 maxcols 1 spacing 1
set xrange [0.0001:]
set grid
set logscale y
set logscale x

f(x) = 0.0000000004*x**2
g(x) = 0.00000001*x**4

plot "builds/area.dat" using 2:3 title "trapezis" with lines, \
     f(x) title "Model Trapezis" with dots, \
     g(x) title "Model Simpson" with dots, \
     "builds/area.dat" using 2:4 title "simpson" with lines
