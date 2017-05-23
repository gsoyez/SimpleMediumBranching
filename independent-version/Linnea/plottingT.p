set style histogram rowstacked gap 0
set style fill solid 0.5 border lt -1
cutoff=0.01
alpha=33.2553/(2*sqrt(0.5)) #Takes x dependence into account, using x=0.5

f(x)=alpha*exp(-alpha*x);
plot "./histogramT.dat" using (0.5*($1+$2)):3 smooth freq with boxes lc rgb"green" t"simulation",\
f(x) t"theory" lt 1 lc 1 lw 2 with lines


