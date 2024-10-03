set output 'figE1.png'
set term png
set xzeroaxis
set yzeroaxis

set grid

f(x) = exp(-x**2/(2*sqrt(2)**2))/(sqrt(2)*sqrt(2*pi))

plot "builds/resE1.dat" index 0 using 1:2:3 title "Histograma" with yerrorbars, \
     f(x) with lines
