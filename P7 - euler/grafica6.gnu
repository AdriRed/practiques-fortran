set key bottom right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

set terminal png 
set output 'P7-23-24-fig6.png'

set xlabel "temps (s)"
set ylabel "energia total (J)"

plot "builds/P7-23-24-res.dat" index 6 using 1:2 with lines title "300", \
    "builds/P7-23-24-res.dat" index 7 using 1:2 with lines title "550", \
    "builds/P7-23-24-res.dat" index 8 using 1:2 with lines title "1000", \
    "builds/P7-23-24-res.dat" index 9 using 1:2 with lines title "20000", \
