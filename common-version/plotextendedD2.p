reset
set terminal pdfcairo dashed enhanced font "Helvetica,6.8"
set title ''

set ylabel 'x* D^2(x,x)'

set output "D2.pdf"

set logscale x
set logscale y

set yrange [1e-4:*]

res_simple="./res/extended-tmax1.0-eps1e-9-xmin1e-4-simple.res"
res="./res/extended-tmax1.0-eps1e-9-xmin1e-4.res"
exact="./res/extended-tmax1.0-exact-xmin1e-4-simple.res"


plot res_simple i 10 using (exp(-$2)):( exp(-$2)*exp(-$2)/(exp(-$1)-exp(-$3)) *$4) with lines lc rgb"blue" dt 1 lw 1  title 'simple',\
res i 10 using (exp(-$2)):( exp(-$2)*exp(-$2)/(exp(-$1)-exp(-$3))*$4) with lines lc rgb"green" dt 2 lw 1  title 'full',\
exact i 10 using (exp(-$2)):((exp(-$2))*$4) with lines dt 3 lc rgb"red" title 'exact'
