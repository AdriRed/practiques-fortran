set key top left
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

set term png
set output "Exa-jan-24-fig1.png"

set xlabel "z"
set ylabel "\phi"

plot "builds/Exa-jan-24-res1.dat" index 0 using 2:3 with lines title "z_0 = 0.3", \
    "builds/Exa-jan-24-res1.dat" index 2 using 2:3 with lines title "z_0 = 0.9"