reset
set terminal pdfcairo enhanced font "Helvetica,13" size 15 cm,9 cm
set title '{Energy loss for different times'
set output "energyloss.pdf"
set key top left


set ylabel 'Loss, average and fluctuation'
set xlabel '{/Symbol \t}'
set xrange [0.1:1.5]

myfile = "./EnergyLoss/energyloss_maxt1.5_xmine-4_itere3epsilone-8.dat"

plot myfile using 1:2 with lines lw 2 lc rgb"blue" title 'th average simple',\
myfile using 1:(sqrt($3)) with lines lw 2 lc rgb"black" title 'th std dev simple',\
myfile using 1:4 with lines lw 3 dt 4 lc rgb"red" title 'num average simple',\
myfile using 1:(sqrt($5)) with lines lw 3 dt 2 lc rgb"red" title 'num std dev simple',\
myfile using 1:6 with lines lw 2 lc rgb"#008000" title 'num average full',\
myfile using 1:(sqrt($7)) with lines lw 2 title 'num std dev full'
