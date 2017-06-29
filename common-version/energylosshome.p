reset
set terminal pdfcairo enhanced dashed font "Helvetica,9" size 7.0, 4.5
set title ''#'{Energy loss for different times'
set output "energyloss.pdf"
set key top left maxrows 3 width -15 at screen .52, screen .46


set ylabel 'Loss, average and fluctuation'
set xlabel '{/Symbol \t}'
set xrange [0.1:1.5]

myfile = "./EnergyLoss/energyloss_maxt1.5_xmine-4_itere3epsilone-8.dat"

plot myfile using 1:2 with lines lt 1 lw 3 lc rgb"red"  title 'Simple (analyt), average',\
myfile using 1:4 with lines lc rgb"red" lt 2 lw 6 title 'Simple (num), average',\
myfile using 1:6 with lines lc rgb"red" lt 3 lw 4  title 'Full (num), average',\
myfile using 1:(sqrt($3)) with lines lt 1 lw 3 lc rgb"blue"  title 'fluctuation',\
myfile using 1:(sqrt($5)) with lines lc rgb"blue" lt 2 lw 6 title 'fluctuation',\
myfile using 1:(sqrt($7)) with lines lc rgb"blue" lt 3 lw 4 title 'fluctuation'
