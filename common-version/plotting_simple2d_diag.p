set style histogram rowstacked gap 0
set style fill solid 0.5 border lt -1
set samples 10000
set logscale x
set xrange [*:1]
set ylabel "corr"
set xlabel "x"
tau=0.1
cutoff=1e-5;
f(x)=( x>cutoff && x<1 ) ? 1/(2*pi)* sqrt(x*x)* 1/sqrt(x*x*(1-x-x)) * (exp((-pi*tau*tau)/(1-x-x)) - exp((-4*pi*tau*tau)/(1-x-x))):0/0; #Nothing if sum above 1
plot "./Results/histogram2d_diag_simple.dat" using (0.5*($1+$2)):3 smooth freq with boxes lc rgb"blue" t"simulation",\
f(x) t"theory" lt 1 lc rgb"red" lw 2 with lines


