reset
set terminal pdfcairo dashed enhanced font "Helvetica,13" #size 10.0,10.0
set title ''
set output "plotD.pdf"
set ylabel '{/Symbol \326}(x) * D(x)' offset 4,0,0
set xlabel 'x'
set key left top maxrows 5

set samples 10000
set logscale x
set logscale y
set xrange [0.01:1]
set yrange [0.001:*]

set lmargin at screen 0.15
set rmargin at screen 0.96

set style line 1 lw 2 lc rgb "red" dt 1
set style line 2 lw 2 lc rgb "black" dt 4
set style line 3 lw 2 lc rgb "blue" dt 1
set style line 4 lw 2 lc rgb "#008000" dt 1
set style line 5 lw 2 lc rgb "#800080" dt 4
set style increment user

D(x,tau)=sqrt(x)* (tau)/(x**0.5 * (1-x)**1.5) * exp((-pi*tau*tau)/(1-x));
D2(x,y,tau)= 1/(2*pi)* sqrt(x*y)* 1/sqrt(x*y*(1-x-y)) * (exp((-pi*tau*tau)/(1-x-y)) - exp((-4*pi*tau*tau)/(1-x-y)));

plot for [tau in ('0.1 0.3 0.5 1 1.2')] D(x,tau+0.0) with lines title '{/Symbol \t}='.tau

set output "plotD2.pdf"
set ylabel 'x* D^2(x,x)'
plot for [tau in ('0.1 0.3 0.6 1 1.2')] D2(x,x,tau+0.0) with lines lw 5 title '{/Symbol \t}='.tau
