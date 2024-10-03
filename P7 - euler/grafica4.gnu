set key bottom right
set key box vertical width 1 height 1 maxcols 1 spacing 1
set grid

set terminal png 
set output 'P7-23-24-fig4.png'

set xlabel "temps (s)"
set ylabel "energia (J)"

plot "builds/P7-23-24-res.dat" index 4 using 1:2 with lines title "Cinetica Euler", \
"builds/P7-23-24-res.dat" index 4 using 1:4 with lines title "Cinetica Euler 2", \
"builds/P7-23-24-res.dat" index 4 using 1:($2+$3) with lines title "Total Euler", \
"builds/P7-23-24-res.dat" index 4 using 1:($4+$5) with lines title "Total Euler 2"
