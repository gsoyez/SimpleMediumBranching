reset
set terminal pdfcairo enhanced font "Helvetica,13" size 15 cm,9 cm
set title '{Energy loss, ratio to theoretical xmin = 0'
set output "energyloss_xmins.pdf"
set key top right


set ylabel 'Simple, num/th'
set xlabel '{/Symbol \t}'
set xrange [0.1:0.7]

myfile = "./EnergyLoss/energyloss_xmins2.dat"

plot 1.0 with lines lt 0 lw 3 lc rgb"black" notitle,\
myfile i 0 using 1:($4/$2) with lines lw 2 title 'xmin = e-1',\
myfile i 1 using 1:($4/$2) with lines lw 2 title 'xmin = e-2',\
myfile i 2 using 1:($4/$2) with lines lw 2 title 'xmin = e-3',\
myfile i 3 using 1:($4/$2) with lines lw 2 lc rgb"red" title 'xmin = e-4',\
myfile i 4 using 1:($4/$2) with lines lw 2 title 'xmin = e-5'


