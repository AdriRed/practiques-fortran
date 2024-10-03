unset key
# set grid

# set xrange [0:100]
set term png
set output "P9-23-24-fig3.png"

set xrange [-0.5:19]
set yrange [-0.5:31.5]
set cbrange [2:35]

set xlabel "x"
set ylabel "y"

plot "builds/data.dat" index 6 using 1:2:3 with image
