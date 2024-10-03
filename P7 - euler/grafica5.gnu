set key bottom right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

set terminal png 
set output 'P7-23-24-fig5.png'

set xlabel "angle (rad)"
set ylabel "velocitat angular (rad)"

plot "builds/P7-23-24-res.dat" index 5 using 2:3 with lines title "Euler 2+", \
"builds/P7-23-24-res.dat" index 5 using 4:5 with lines title "Euler 2-", \
