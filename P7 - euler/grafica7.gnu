unset key

set terminal gif animate delay 0.00001
set output 'P7-23-24-figextra.gif'
set xlabel "x (m/L)"
set ylabel "y (m/L)"
set grid

set xrange [-1:1]
set yrange [-1:1]

do for [i=10:1010] {
    plot 'builds/P7-23-24-res.dat' index (i) using (cos($2+pi/2)):(-sin($2+pi/2)) with circles
}