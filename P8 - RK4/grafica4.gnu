set key top left
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set xrange [-7:1]
set term png
set output "P8-23-24-fig4.png"

set xlabel "x (\angstrom)"
set ylabel "\phi (\angstrom^{-1/2})"

plot \
    "builds/P8-23-24-res.dat" index 10 using 1:2 with lines title "\beta = 0", \
    "builds/P8-23-24-res.dat" index 11 using 1:2 with lines title "\beta = 5", \
    "builds/P8-23-24-res.dat" index 12 using 1:2 with lines title "\beta = 15", \
