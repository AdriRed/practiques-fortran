set key top right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid
set format x "%.tx10^{%L}"

set term png
set output "P6-23-24-fig3.png"

set xlabel "Iteracions"
set ylabel "Error"

plot "builds/P6-23-24-res.dat" index 2 using 1:3 with linespoints title "Error estimat", \
    "builds/P6-23-24-res.dat" index 2 using 1:(abs($2-pi)) with linespoints title "Error real"