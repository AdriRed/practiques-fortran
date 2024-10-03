set key bottom left
set key box vertical width 1 height 1 maxcols 1 spacing 1

set logscale xy

set xrange [1.6e3:1.3e5]


I_1 = pi**3/2.
I_2 = 6*pi**3 - 905.*pi/144.

plot "builds/data.dat" index 0 using 1:(abs($2-I_1)) with lines title "Error real I_1", \
    "builds/data.dat" index 0 using 1:3 title "Error estimat I_1", \
    "builds/data.dat" index 0 using 1:(abs($4-I_2)) with lines title "Error real I_2", \
    "builds/data.dat" index 0 using 1:5 title "Error estimat I_2"