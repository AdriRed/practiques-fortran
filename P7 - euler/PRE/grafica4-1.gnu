set key bottom right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set terminal png 
# set output 'grafica1.png'

set xlabel "temps (s)"
set ylabel "energia potencial (J)"

plot "builds/P7-23-24-res.dat" index 4 using 1:3 with points title "Potencial Euler (1s)", \
"builds/P7-23-24-res.dat" index 4 using 1:5 with points title "Potencial Adam-Bashford (1s)"
