set key center right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set xrange [-7:1]
set term png
set output "P8-23-24-fig2.png"

set xlabel "Iteracions"
set ylabel "Energia (eV)"

plot \
    "builds/P8-23-24-res.dat" index 4 using 1:2 with linespoints title "E_1", \
    "builds/P8-23-24-res.dat" index 6 using 1:2 with linespoints title "E_2", \
    "builds/P8-23-24-res.dat" index 8 using 1:2 with linespoints title "E_3"
