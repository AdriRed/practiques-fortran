set key bottom right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set terminal png 
# set output 'grafica1.png'

# set xrange [-1:23]
# set yrange [-0.015:0.015]

set xlabel "angle (rad)"
set ylabel "velocitat angular (rad)"

plot "builds/P7-23-24-res.dat" index 6 using 2:3 with lines title "Adams-Bashford+", \
"builds/P7-23-24-res.dat" index 6 using 4:5 with lines title "Adams-Bashford-", \
