set xlabel 'E'
set output 'P3-23-24-fig3.png'
set term png
set xzeroaxis
set yzeroaxis
set key left top
set key box vertical width 1 height 1 maxcols 1 spacing 1
set yrange [-0.5:1]
set xrange [0:2]
set grid

plot "builds/P3-23-24-res3-n25.dat" using 1:3 with lines title "Approx 25", \
     "builds/P3-23-24-res3-n230.dat" using 1:3 with lines title "Approx 230", \
     "builds/P3-23-24-res3-n230.dat" using 1:4 pointsize 0.5 title "dF"

