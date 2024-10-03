set key top left
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

# set xrange [-7:1]
set term png
set output "P8-23-24-fig3.png"

set xlabel "x (\angstrom)"
set ylabel "\phi (\angstrom^{-1/2})"

plot \
    "builds/P8-23-24-res.dat" index 5 using 1:2 with lines title "\phi_1", \
    "builds/P8-23-24-res.dat" index 7 using 1:2 with lines title "\phi_2", \
    "builds/P8-23-24-res.dat" index 9 using 1:2 with lines title "\phi_3", \
