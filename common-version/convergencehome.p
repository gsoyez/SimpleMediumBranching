reset
set terminal pdfcairo enhanced font "Helvetica,7"
set title '{/Symbol \326}(x) * D(x) at {/Symbol \t}=1'
set output "convergence.pdf"
set logscale x
set logscale y; set ytics (0.8,0.95,1,1.05,""1.2,""1.4,""1.6,""1.8,2)
set yrange [0.8:2]
set bars small

set ylabel 'Simple, ratio numerical/analytic' offset 2,0,0
set xlabel 'x'

res_5='<bash -c "paste <(grep -v \"^#\" ./res/extended--tmax1.0-eps1e-5-xmin1e-4.res) <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-9-xmin1e-4.res)"'
res_6='<bash -c "paste <(grep -v \"^#\" ./res/extended--tmax1.0-eps1e-6-xmin1e-4.res) <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-9-xmin1e-4.res)"'
res_7='<bash -c "paste <(grep -v \"^#\" ./res/extended--tmax1.0-eps1e-7-xmin1e-4.res) <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-9-xmin1e-4.res)"'
res_8='<bash -c "paste <(grep -v \"^#\" ./res/extended--tmax1.0-eps1e-8-xmin1e-4.res) <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-9-xmin1e-4.res)"'
res_9='<bash -c "paste <(grep -v \"^#\" ./res/extended--tmax1.0-eps1e-9-xmin1e-4.res) <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-9-xmin1e-4.res)"'

res_5_simple='<bash -c "paste <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-5-xmin1e-4-simple.res) <(grep -v \"^#\" ./res/extended-tmax1.0-exact-xmin1e-4-simple.res)"'
res_6_simple='<bash -c "paste <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-6-xmin1e-4-simple.res) <(grep -v \"^#\" ./res/extended-tmax1.0-exact-xmin1e-4-simple.res)"'
res_7_simple='<bash -c "paste <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-7-xmin1e-4-simple.res) <(grep -v \"^#\" ./res/extended-tmax1.0-exact-xmin1e-4-simple.res)"'
res_8_simple='<bash -c "paste <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-8-xmin1e-4-simple.res) <(grep -v \"^#\" ./res/extended-tmax1.0-exact-xmin1e-4-simple.res)"'
res_9_simple='<bash -c "paste <(grep -v \"^#\" ./res/extended-tmax1.0-eps1e-9-xmin1e-4-simple.res) <(grep -v \"^#\" ./res/extended-tmax1.0-exact-xmin1e-4-simple.res)"'

exact="./res/extended-tmax1.0-exact-xmin1e-4-simple.res"

plot res_5_simple i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) with yerr ps 0 lw 3 title '{/Symbol \e}=10^{-5}',\
res_6_simple i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) with yerr ps 0 lw 3 title '{/Symbol \e}=10^{-6}',\
res_7_simple i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) with yerr ps 0 lw 3 title '{/Symbol \e}=10^{-7}',\
res_8_simple i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) with yerr ps 0 lw 3 title '{/Symbol \e}=10^{-8}',\
res_9_simple i 9 using (exp(-$2)):(sqrt(exp(-$2))*$4)/(sqrt(exp(-$2))*$9):(sqrt(exp(-$2))*$5)/(sqrt(exp(-$2))*$9) with yerr ps 0 lw 3 title '{/Symbol \e}=10^{-9}'


