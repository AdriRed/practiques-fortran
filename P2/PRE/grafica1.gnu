set xlabel 't (s)'
set ylabel 'x (cm)'
set output 'P2-23-24-fig1.png'
set term png
set xzeroaxis
set yzeroaxis
set key right bottom

plot "builds/P2-23-24-res1.dat" using 1:3 title "P2", "builds/P2-23-24-res1.dat" using 1:4  title "P3"
    
