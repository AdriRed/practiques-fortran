set xlabel 't (years)'
set ylabel 'Area (AU^2)'
set output 'P3-23-24-fig3.png'
set term png
set xzeroaxis
set yzeroaxis
unset key
# set yrange [-0.5:1]
# set xrange [0:2*pi]
set grid

plot "builds/P3-23-24-res.dat" index 3 using 1:2
