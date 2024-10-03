set xlabel 'N'
set ylabel 'S^8_N'
set output 'P1-23-24-fig1.png'
set term png
set xzeroaxis
set yzeroaxis
set key right bottom

plot "builds/P1-23-24-res1.dat" using 1:3  title "Model", "builds/P1-23-24-res1.dat" using 1:2 title "S^M_N" pointsize 1
    
