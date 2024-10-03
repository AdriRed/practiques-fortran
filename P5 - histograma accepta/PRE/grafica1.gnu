set ylabel 'x'
set xlabel 'distribucio'
# set key left top
# set key box vertical width 1 height 1 maxcols 1 spacing 1
set xrange [-pi:pi]
# set style data histogram

# set output 'P5-23-24-b-fig1.png'
# set terminal png

p(x) = 125.*exp(pi)*(x*sin(x))**2*exp(-abs(x))/(4.*(68.*exp(pi)-70.*pi-25.*pi**2-68.))

plot p(x), \
    "builds/P5-23-24.dat" index 1 using 1:2:3 with boxes title "Histograma"