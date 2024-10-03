set key center right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set xrange [0:400]
set term png
set output "P9-23-24-fig3.png"

set xlabel "Iteracions"
set ylabel "Temperatura (ºC)"

set title "T_{ini} = 1040ºC"

plot \
    "builds/P9-23-24-res.dat" index 4 using 1:2 with lines title "Gauss Seidel", \
    "builds/P9-23-24-res.dat" index 5 using 1:2 with lines title "Sobrerelaxació", 