reset
set terminal pdfcairo dashed enhanced font "Helvetica,6.8" size 5.0,7.0
set title ''

set ylabel 'x* D^2(x,x)'

set output "D2.pdf"

set key left top maxrows 3


set multiplot

set lmargin at screen 0.15
set rmargin at screen 0.96

set bmargin at screen 0.55
set tmargin at screen 0.96


set logscale x
set logscale y

res_simple="./res/extended-tmax1.0-eps1e-9-xmin1e-4-simple.res"
res="./res/extended-tmax1.0-eps1e-9-xmin1e-4.res"
exact="./res/extended-tmax1.0-exact-xmin1e-4-simple.res"

set yrange [0.001:*]

plot res_simple i 0 using (exp(-$2)):(sqrt(exp(-$2))*$4) with lines lc rgb"red" lt 2 lw 5  title 'simple {/Symbol \t}=0.1',\
res i 0 using (exp(-$2)):(sqrt(exp(-$2))*$4) with lines lc rgb"red" lt 2 lw 3  title 'full {/Symbol \t}=0.1',\
exact i 0 using (exp(-$2)):(sqrt(exp(-$2))*$4) with lines lt 1 lc rgb"red" title 'exact {/Symbol \t}=0.1',\
res_simple i 4 using (exp(-$2)):(sqrt(exp(-$2))*$4) with lines lc rgb"blue" lt 3 lw 5  title 'simple {/Symbol \t}=0.4',\
res i 4 using (exp(-$2)):(sqrt(exp(-$2))*$4):(sqrt(exp(-$2))*$5) with lines lc rgb"blue" lt 3 lw 3   title 'full {/Symbol \t}=0.4',\
exact i 4 using (exp(-$2)):(sqrt(exp(-$2))*$4) with lines lt 1 lc rgb"blue" title 'exact {/Symbol \t}=0.4',\
res_simple i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4):(sqrt(exp(-$2))*$5) with lines lc rgb"#008000" lt 4 lw 5  title 'simple {/Symbol \t}=1',\
res i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4)with lines lc rgb"#008000" lt 4 lw 3   title 'full {/Symbol \t}=1',\
exact i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4) with lines lt 1 lc rgb"#008000" title 'exact {/Symbol \t}=1'





combined='<bash -c "paste <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-9-xmin1e-4-simple.res) <(grep -v \"^#\" ./res/extended-tmax1.0-exact-xmin1e-4-simple.res)"'


unset ylabel
set ylabel "Relative error"

set bmargin at screen 0.40
set tmargin at screen 0.50
unset logscale y
set yrange [0.95:1.05]

plot 1.0 lc rgb"black" lt 0 lw 3 notitle,\
combined i 0 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) t 'simple {/Symbol \t}=0.1' w yerr lc rgb"red" lw 1 pt 7 ps 0 lt 1


set bmargin at screen 0.25
set tmargin at screen 0.35
unset logscale y
set yrange [0.95:1.05]

plot 1.0 lc rgb"black" lt 0 lw 3 notitle,\
combined i 4 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) t 'simple {/Symbol \t}=0.4' w yerr lc rgb"red" lw 1 pt 7 ps 0 lt 1


set xlabel "x"

set bmargin at screen 0.10
set tmargin at screen 0.20
unset logscale y
set yrange [0.95:1.05]

plot 1.0 lc rgb"black" lt 0 lw 3 notitle,\
combined i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) t 'simple {/Symbol \t}=1' w yerr lc rgb"red" lw 1 pt 7 ps 0 lt 1







unset multiplot
     
