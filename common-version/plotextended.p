reset
set terminal pdfcairo enhanced font "Helvetica,11" size 5.0,6.0
set title ''

set ylabel '{/Symbol \326}(x) * D(x)'

set output "times.pdf"
set key right top outside


set multiplot

set lmargin at screen 0.15
set rmargin at screen 0.8

set bmargin at screen 0.55
set tmargin at screen 0.96


set logscale x
set logscale y
res_simple="./res/extended-tmax1.0-eps1e-9-xmin1e-4-simple.res"
res="./res/extended-tmax1.0-eps1e-9-xmin1e-4.res"
exact="./res/extended-tmax1.0-exact-xmin1e-4-simple.res"

set yrange [0.001:*]

plot res_simple i 0 using (exp(-$2)):(sqrt(exp(-$2))*$4):(sqrt(exp(-$2))*$5) with yerr lc rgb"#8B0000" pt 0 ps 0  title 'simple {/Symbol \t}=0.1',\
res i 0 using (exp(-$2)):(sqrt(exp(-$2))*$4):(sqrt(exp(-$2))*$5) with yerr lc rgb"red" pt 0 ps 0  title 'full {/Symbol \t}=0.1',\
exact i 0 using (exp(-$2)):(sqrt(exp(-$2))*$4) with lines lt 1 lc rgb"orange" title 'exact {/Symbol \t}=0.1',\
res_simple i 4 using (exp(-$2)):(sqrt(exp(-$2))*$4):(sqrt(exp(-$2))*$5) with yerr lc rgb"black"  pt 0 ps 0  title 'simple {/Symbol \t}=0.4',\
res i 4 using (exp(-$2)):(sqrt(exp(-$2))*$4):(sqrt(exp(-$2))*$5) with yerr lc rgb"#4169E1" pt 0 ps 0  title 'full {/Symbol \t}=0.4',\
exact i 4 using (exp(-$2)):(sqrt(exp(-$2))*$4) with lines lt 1 lc rgb"#00FFFF" title 'exact {/Symbol \t}=0.4',\
res_simple i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4):(sqrt(exp(-$2))*$5) with yerr lc rgb"#006400"  pt 0 ps 0  title 'simple {/Symbol \t}=1',\
res i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4):(sqrt(exp(-$2))*$5) with yerr lc rgb"#6B8E23" pt 0 ps 0  title 'full {/Symbol \t}=1',\
exact i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4) with lines lt 1 lc rgb"#32CD32" title 'exact {/Symbol \t}=1'











combined='<bash -c "paste <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-9-xmin1e-4-simple.res) <(grep -v \"^#\" ./res/extended-tmax1.0-exact-xmin1e-4-simple.res)"'


unset ylabel
set ylabel "Relative error"

set bmargin at screen 0.40
set tmargin at screen 0.50
unset logscale y
set yrange [0.9:1.1]

plot 1.0 lc rgb"black" lt 0 lw 3,\
combined i 0 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) t 'simple {/Symbol \t}=0.1' w yerr lc rgb"red" lw 1 pt 7 ps 0


set bmargin at screen 0.25
set tmargin at screen 0.35
unset logscale y
set yrange [0.9:1.1]

plot 1.0 lc rgb"black" lt 0 lw 3,\
combined i 4 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) t 'simple {/Symbol \t}=0.4' w yerr lc rgb"red" lw 1 pt 7 ps 0


set xlabel "x"

set bmargin at screen 0.10
set tmargin at screen 0.20
unset logscale y
set yrange [0.9:1.1]

plot 1.0 lc rgb"black" lt 0 lw 3,\
combined i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) t 'simple {/Symbol \t}=1' w yerr lc rgb"red" lw 1 pt 7 ps 0







unset multiplot
     
