set xlabel 'x_1 (cm)'
set ylabel 'x (cm)'
set output 'P2-23-24-fig2-c.png'
set term png
set xzeroaxis
set yzeroaxis
set key right bottom

set key box vertical width 1 height 1 maxcols 1 spacing 1

plot "builds/P2-23-24-res1-c.dat" using 2:4 title "P2", \
     "builds/P2-23-24-res1-c.dat" using 2:6 title "P3"
    
