unset key
# set grid

# set xrange [0:100]
stats 'builds/P10-23-24-res.dat'
set term gif animate delay 0.01
set output "animation.gif"

set xlabel "x"
set ylabel "T"

do for [i=3:int(STATS_blocks-1)]{
    plot 'builds/P10-23-24-res.dat' index i using 1:2 with linespoints
}