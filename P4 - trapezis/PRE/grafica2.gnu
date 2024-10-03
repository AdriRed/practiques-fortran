set ylabel 'Divergencia massa'
set xlabel 'h'
# set xzeroaxis
set key left top
set key box vertical width 1 height 1 maxcols 1 spacing 1
# set yrange [-0.5:1]
set grid

plot "builds/massa.dat" using 2:3 title "trapezis" with lines, \
     "builds/massa.dat" using 2:4 title "simpson" with lines
