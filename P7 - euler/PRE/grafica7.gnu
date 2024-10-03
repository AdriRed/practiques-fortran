set key bottom right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set terminal png 
# set output 'grafica1.png'

# set xrange [-1:23]
# set yrange [-0.015:0.015]

set xlabel "angle (rad)"
set ylabel "velocitat angular (rad)"

plot "builds/P7-23-24-res.dat" index 7 using 1:3 with lines title "300", \
    "builds/P7-23-24-res.dat" index 8 using 1:3 with lines title "1000", \
    "builds/P7-23-24-res.dat" index 9 using 1:3 with lines title "2200", \
    "builds/P7-23-24-res.dat" index 10 using 1:3 with lines title "14500", \
