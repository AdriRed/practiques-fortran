set key top right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid
set format x "%.tx10^{%L}"

set term png
set output "P6-23-24-fig4.png"

set xlabel "Iteracions"

f(x) = pi

plot "builds/P6-23-24-res.dat" index 2 using 1:2:3 with errorbars title "Calcul de pi", \
    f(x) title "pi"