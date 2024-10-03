set key center right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

set term png
set output "P6-23-24-fig1.png"

set xlabel "Iteracions"

plot "builds/P6-23-24-res.dat" index 0 using 1:2:3 with errorbars title "I_{up}", \
    "builds/P6-23-24-res.dat" index 0 using 1:4:5 with errorbars title "I_{down}"
