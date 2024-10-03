set key bottom left
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

set terminal png 
set output 'P7-23-24-fig1.png'


set xlabel "t (s)"
set ylabel "velocitat angular (rad)"

plot "builds/P7-23-24-res.dat" index 0 using 1:3 with lines title "Euler", \
"builds/P7-23-24-res.dat" index 1 using 1:3 with lines title "Euler 2", \
