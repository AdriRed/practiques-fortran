set key bottom right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set terminal png 
# set output 'grafica1.png'

# set xrange [-1:23]
# set yrange [-0.015:0.015]

set xlabel "temps (s)"
set ylabel "energia total (J)"

plot "builds/P7-23-24-res.dat" index 4 using 1:($2+$3) with lines title "Potencial Euler", \
"builds/P7-23-24-res.dat" index 4 using 1:($4+$5) with lines title "Potencial Adam-Bashford", \
