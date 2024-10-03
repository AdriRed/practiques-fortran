unset key
# set grid

# set xrange [0:100]
stats 'builds/animation.dat'
set term gif animate delay 0.01
set output "P9-23-24-fig3.gif"

set xrange [-0.5:19]
set yrange [-0.5:31.5]
set cbrange [2:35]

set xlabel "x"
set ylabel "y"

do for [i=0:int(STATS_blocks-1)]{
    plot "builds/animation.dat" index i using 1:2:3 with image
}
