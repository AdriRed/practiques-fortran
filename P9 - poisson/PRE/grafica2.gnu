set key center right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set xrange [0:100]
set term png
set output "P9-23-24-fig2.png"

set xlabel "Iteracions"
set ylabel "Temperatura (ºC)"

set title "Sobrerelaxacio Successiva"

plot \
    "builds/data.dat" index 1 using 1:2 with lines title "T_1", \
    "builds/data.dat" index 3 using 1:2 with lines title "T_2", \
    "builds/data.dat" index 5 using 1:2 with lines title "T_3"
