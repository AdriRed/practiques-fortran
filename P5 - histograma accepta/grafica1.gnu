set key center top
set key box vertical width 1 height 1 maxcols 1 spacing 1
# set xrange [-pi:pi]
# set style data histogram

# set output 'P5-23-24-b-fig1.png'
# set terminal png

set xlabel "x (nm)"
set ylabel 'Probabilitat'

L = 4
p(x) = sin(x/L)**2/(L*pi)

plot p(x), \
    "builds/P5-23-24-res.dat" index 0 using 1:2:3 with yerrorbars  title "Histograma"