set xlabel 'x (AU)'
set ylabel 'y (AU)'
set output 'P3-23-24-fig2.png'
set term png
set xzeroaxis
set yzeroaxis
set key left top
set key box vertical width 1 height 1 maxcols 1 spacing 1
# set yrange [-0.5:1]
# set xrange [0:2*pi]
set grid

plot "builds/P3-23-24-res.dat" index 2 using 3:4 title "Halley" \
