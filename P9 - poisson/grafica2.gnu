unset key

set term png
set output "P9-23-24-fig4.png"

set xrange [-0.5:33.5]
set yrange [-0.5:45.5]

set xlabel "x (cm)"
set ylabel "y (cm)"
set clabel "T (ºC)"

plot "builds/P9-23-24-res.dat" index 6 using 1:2:3 with image
