reset
set terminal pdfcairo dashed enhanced font "Helvetica,7" size 5.0,7.0
set title ''

set ylabel 'x D^{(2)}(x,x)' offset 4,0,0

set output "D2.pdf"

set key maxrows 3 width -21 top right at screen .87, screen .95 

set multiplot

set lmargin at screen 0.10
set rmargin at screen 0.96

set bmargin at screen 0.55
set tmargin at screen 0.96


set logscale x
set logscale y

res_simple="./res/extended-tmax1.0-eps1e-9-xmin1e-4-simple.res"
res="./res/extended-tmax1.0-eps1e-9-xmin1e-4.res"
exact="./res/extended-tmax1.0-exact-xmin1e-4-simple.res"

set yrange [0.001:*]
set yrange [0.0001:1]

D2(x,y,tau)= 1/(2*pi)* sqrt(x*y)* 1/sqrt(x*y*(1-x-y)) * (exp((-pi*tau*tau)/(1-x-y)) - exp((-4*pi*tau*tau)/(1-x-y)));


plot res_simple i 10 using (exp(-$2)):(D2(exp(-$2),exp(-$2),0.1)) with lines lt 1 lw 3 lc rgb"red" title 'simple, analytic:   {/Symbol \t}=0.1',\
res_simple i 10 using (exp(-$2)):((exp(-$2))*$4) with lines lc rgb"red" lt 2 lw 6  title 'simple, numerical:   {/Symbol \t}=0.1',\
res i 10 using (exp(-$2)):((exp(-$2))*$4) with lines lc rgb"red" lt 3 lw 4  title 'full, numerical:   {/Symbol \t}=0.1',\
res_simple i 14 using (exp(-$2)):(D2(exp(-$2),exp(-$2),0.5)) with lines lt 1 lw 3 lc rgb"blue" title '{/Symbol \t}=0.5',\
res_simple i 14 using (exp(-$2)):((exp(-$2))*$4) with lines lc rgb"blue" lt 2 lw 6  title '{/Symbol \t}=0.5',\
res i 14 using (exp(-$2)):((exp(-$2))*$4) with lines lc rgb"blue" lt 3 lw 4  title '{/Symbol \t}=0.5',\
res_simple i 19 using (exp(-$2)):(D2(exp(-$2),exp(-$2),1)) with lines lt 1 lw 3 lc rgb"#008000" title '{/Symbol \t}=1',\
res_simple i 19 using (exp(-$2)):((exp(-$2))*$4):((exp(-$2))*$5) with lines lc rgb"#008000" lt 2 lw 6  title '{/Symbol \t}=1',\
res i 19 using (exp(-$2)):((exp(-$2))*$4)with lines lc rgb"#008000" lt 3 lw 4   title '{/Symbol \t}=1'




unset ylabel
set key default
set key top left

set bmargin at screen 0.40
set tmargin at screen 0.50
unset logscale y
set yrange [0.95:1.05]
plot 1.0 lc rgb"black" lt 0 lw 3 notitle,\
res_simple i 10 using (exp(-$2)):((exp(-$2))*$4)/(D2(exp(-$2),exp(-$2),0.1)):((exp(-$2))*$5)/(D2(exp(-$2),exp(-$2),0.1)) t 'simple {/Symbol \t}=0.1' w yerr lc rgb"red" lw 3 pt 7 ps 0 lt 1

set ylabel "Simple, ratio numerical/analytic" offset 2,0,0


set bmargin at screen 0.25
set tmargin at screen 0.35
unset logscale y
#set yrange [0.95:1.05]

plot 1.0 lc rgb"black" lt 0 lw 3 notitle,\
res_simple i 14 using (exp(-$2)):((exp(-$2))*$4)/(D2(exp(-$2),exp(-$2),0.5)):((exp(-$2))*$5)/(D2(exp(-$2),exp(-$2),0.5)) t 'simple {/Symbol \t}=0.5' w yerr lc rgb"red" lw 3 pt 7 ps 0 lt 1

unset ylabel
set xlabel "x"

set bmargin at screen 0.10
set tmargin at screen 0.20
unset logscale y
#set yrange [0.95:1.05]

plot 1.0 lc rgb"black" lt 0 lw 3 notitle,\
res_simple i 19 using (exp(-$2)):((exp(-$2))*$4)/(D2(exp(-$2),exp(-$2),1)):((exp(-$2))*$5)/(D2(exp(-$2),exp(-$2),1)) t 'simple {/Symbol \t}=1' w yerr lc rgb"red" lw 3 pt 7 ps 0 lt 1



unset multiplot
     
