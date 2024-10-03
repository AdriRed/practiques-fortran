set key top right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set terminal png 
# set output 'grafica1.png'

# set xrange [-1:23]
# set yrange [-0.015:0.015]

set xlabel "t (s)"
set ylabel "velocitat angular (rad)"

plot "builds/P7-23-24-res.dat" index 2 using 1:3 with lines title "Euler", \
"builds/P7-23-24-res.dat" index 3 using 1:3 with lines title "Adam-Bashford", \
