set xlabel 'x_3 (s)'
set ylabel 'x (cm)'
set output 'P2-23-24-fig2.png'
set term png
set xzeroaxis
set yzeroaxis
set key right bottom

plot "builds/P2-23-24-res1.dat" using 4:2 title "P1", "builds/P2-23-24-res1.dat" using 4:5  title "P4"
    
