set ylabel 'Divergencia area (10^{12} km^2)'
set xlabel 'h (10^6 km)'
set key left top
set key box vertical width 1 height 1 maxcols 1 spacing 1
set yrange [1e-12:1e2]
set grid
set logscale y
set logscale x

set format xy "10^{%T}"

set output 'P4-23-24-b-fig1.png'
set terminal png

f(x) = 0.03*x**2
g(x) = 0.0000001*x**4

plot "builds/P4-23-24-res.dat" using 1:4 index 0 title "Error Trapezis" with points, \
     f(x) title "Model Trapezis" with lines, \
     "builds/P4-23-24-res.dat" using 1:5 index 0 title "Error Simpson" with points, \
     g(x) title "Model Simpson" with lines
