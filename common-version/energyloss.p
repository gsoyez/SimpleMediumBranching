reset
set terminal pdfcairo enhanced font "Helvetica,11"
set title '{Energy loss for different times'
set output "energyloss.pdf"
set key top left


set ylabel 'Loss, average and fluctuation'
set xlabel '{/Symbol \t}'

myfile = "./EnergyLoss/energyloss_full_e-4.dat"

plot myfile using 1:2 with lines title 'th average',\
myfile using 1:3 with lines title 'num average',\
myfile using 1:(sqrt($4)) with lines title 'th std dev',\
myfile using 1:(sqrt($5)) with lines title 'num std dev'
