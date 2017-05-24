set logscale x
set xlabel 'x'
set ylabel 'Fraction de gluons'
set xrange[0.0000002:1]
set grid
set nokey

plot "output.dat" using 1:2 with lines 

pause -1 

