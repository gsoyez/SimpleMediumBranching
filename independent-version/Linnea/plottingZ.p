set style histogram rowstacked gap 0
set style fill solid 0.5 border lt -1
cutoff=0.01
alpha=33.2553/2 #Here alpha is just to normalize k(z) => indep of x

f(x)=(x>cutoff&&x<1-cutoff) ? (1/(alpha))*((1-x*(1-x))**2.5)/((x*(1-x))**1.5):0/0;
plot "./histogramZ.dat" using (0.5*($1+$2)):3 smooth freq with boxes lc rgb"green" t"simulation",\
f(x) t"theory" lt 1 lc 1 lw 2 with lines


