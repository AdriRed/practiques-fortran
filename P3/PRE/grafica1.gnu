set xlabel 'E'
set output 'P3-23-24-fig1.png'
set term png
set xzeroaxis
set yzeroaxis
set key left top
set key box vertical width 1 height 1 maxcols 1 spacing 1
set yrange [-0.5:1]
set xrange [0:2.5]
set grid

plot "builds/P3-23-24-res.dat" index 0 using 1:2 title "f(E)" with lines, \
     "builds/P3-23-24-res.dat" index 0 using 1:3 title "df(E)" with lines
