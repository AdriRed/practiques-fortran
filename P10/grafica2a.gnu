unset key
# set grid

# set xrange [0:100]
stats 'builds/P10-23-24-res-2a.dat'

set term png
set output "grafica2a.png"

set xlabel "t (s)"
set ylabel "T - T_{amb} (ºC)"

plot 'builds/P10-23-24-res-2a.dat' using 1:2 every ::3::3