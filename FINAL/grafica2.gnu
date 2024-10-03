unset key
set grid

set term png
set output "Exa-jan-24-fig2.png"
set xrange[0.5: 10.5]
set xlabel "N (*10^4)"
set ylabel "Integral"

plot "builds/Exa-jan-24-res1.dat" index 5 using 1:2:3 with yerrorbars 