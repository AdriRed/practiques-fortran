set xlabel 'N'
set ylabel 'Fermi energy (eV)'
set output 'P1-23-24-fig1.png'
set term png
set xzeroaxis
set yzeroaxis

plot "builds/P1-23-24-res1.dat" using 1:2 title "Fermi Energy" pointsize 2
