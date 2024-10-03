set key top right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

set terminal png 
set output 'P7-23-24-fig2.png'

set xlabel "t (s)"
set ylabel "angle (rad)"

plot "builds/P7-23-24-res.dat" index 2 using 1:2 with lines title "Euler", \
"builds/P7-23-24-res.dat" index 3 using 1:2 with lines title "Euler 2", \
