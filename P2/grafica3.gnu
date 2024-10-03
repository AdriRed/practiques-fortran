set ylabel 'x_4 (cm)'
set xlabel 't (s)'
set output 'P2-23-24-fig3-c.png'
set term png
set xzeroaxis
set yzeroaxis
set xrange [0:3]
set key center top
set key box vertical width 1 height 1 maxcols 1

plot "builds/P2-23-24-res2-c.dat" using 1:2 title "Lineal", \
     "builds/P2-23-24-res3-c.dat" using 1:2 title "Ordre 0", \
     "builds/P2-23-24-res1-c.dat" using 1:5 title "Original" 
    
