set ylabel 'Iteracions'
set xlabel 'E_0'
set output 'P3-23-24-fig2.png'
set term png
set xzeroaxis
set yzeroaxis
set nokey
set yrange [0:]
set xrange [0:5.4]
set grid

set style histogram rows
set boxwidth 0.05

plot "builds/P3-23-24-res.dat" index 1 using 1:2 \
     with boxes fill pattern 3
