set key center right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set xrange [0:400]
set term png
set output "P9-23-24-fig1.png"

set xlabel "Iteracions"
set ylabel "Temperatura (ºC)"

set title "Jacobi"

plot \
    "builds/data.dat" index 0 using 1:2 with lines title "T_1", \
    "builds/data.dat" index 2 using 1:2 with lines title "T_2", \
    "builds/data.dat" index 4 using 1:2 with lines title "T_3"
