set xlabel 'x_2 (cm)'
set ylabel 't (s)'
set output 'P2-23-24-fig3.png'
set term png
set xzeroaxis
set yzeroaxis
set key right bottom
set xrange [0:3]

plot "builds/P2-23-24-res2.dat" using 1:2 title "Interpolació", \
    "builds/P2-23-24-res1.dat" using 1:3 title "Original" 
    
