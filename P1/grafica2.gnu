set xlabel 'N'
set ylabel 'S^8_N/S^{asim}_N'
set output 'P1-23-24-fig2.png'
set term png
set xzeroaxis
set yzeroaxis
set key right bottom

plot "builds/P1-23-24-res1.dat" using 1:4 title "Relacion" pointsize 1
    
