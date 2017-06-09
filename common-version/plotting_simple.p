set style histogram rowstacked gap 0
set style fill solid 0.5 border lt -1
set samples 10000
set logscale x
set xrange [*:1]
set ylabel "sqrt(x)*D(x)"
set xlabel "x"
tau=0.1
cutoff=1e-4;
f(x)=(x>cutoff&&x<1) ? sqrt(x)* (tau)/(x**0.5 * (1-x)**1.5) * exp((-pi*tau*tau)/(1-x)):0/0;
plot "./histogram_simple.dat" using (0.5*($1+$2)):3 smooth freq with boxes lc rgb"blue" t"simulation",\
f(x) t"theory" lt 1 lc rgb"red" lw 2 with lines


