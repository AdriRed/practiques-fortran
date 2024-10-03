set xlabel 'N'
set ylabel 'Energy (eV)'
set output 'P1-23-24-fig2.png'
set term png
set xzeroaxis
set yzeroaxis
set xrange [1:20]
set key right bottom

f(x) = 8 - 6/x + 6/x**2

plot f(x) title "Model", "builds/P1-23-24-res1.dat" using 1:3 title "2N/N Fermi Energy" pointsize 2
    
