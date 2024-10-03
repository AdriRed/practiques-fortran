set xlabel 't (s)'
set ylabel 'x (cm)'
set output 'P2-23-24-fig1-c.png'
set term png
set xzeroaxis
set yzeroaxis
set key right bottom
set key box vertical width 1 height 1 maxcols 1 spacing 1
set key outside

plot "builds/P2-23-24-res1-c.dat" using 1:2 title "P1", \
     "builds/P2-23-24-res1-c.dat" using 1:4 title "P3", \
     "builds/P2-23-24-res1-c.dat" using 1:6 title "P5"
    
