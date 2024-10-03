# set key left top
# set key box vertical width 1 height 1 maxcols 1 spacing 1
# set xrange [-pi:pi]
# set style data histogram

# set output 'P5-23-24-b-fig1.png'
# set terminal png

set xlabel 'x (um)'
set ylabel 'Probabilitat'

SIGMA = 4
g(x) =  exp(-x**2./(2.*SIGMA**2))/sqrt(2.*pi*SIGMA)

plot g(x), \
    "builds/P5-23-24-res.dat" index 2 using 1:2:3 with yerrorbars  title "Histograma"