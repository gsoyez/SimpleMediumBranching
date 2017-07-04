reset
set terminal pdfcairo enhanced font "Helvetica,13" size 15 cm,9 cm
set title ''#'{Energy loss for different times'
set output "energyloss.pdf"
set key top left maxrows 3 width -15 at screen .52, screen .46


set ylabel 'Loss, average and fluctuation'
set xlabel '{/Symbol \t}'
set xrange [0.1:2.4]

myfile = "./EnergyLoss/energyloss_maxt2.4_xmine-4_itere3epsilone-8.dat"

plot myfile using 1:2 with lines dt 1 lw 1 lc rgb"red"  title 'Simple (analyt), average',\
myfile using 1:4 with lines lc rgb"red" dt 2 lw 3 title 'Simple (num), average',\
myfile using 1:6 with lines lc rgb"red" dt 3 lw 3  title 'Full (num), average',\
myfile using 1:(sqrt($3)) with lines dt 1 lw 1 lc rgb"blue"  title 'fluctuation',\
myfile using 1:(sqrt($5)) with lines lc rgb"blue" dt 2 lw 3 title 'fluctuation',\
myfile using 1:(sqrt($7)) with lines lc rgb"blue" dt 3 lw 3 title 'fluctuation'
