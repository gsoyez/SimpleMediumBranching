set title "Correlations"
set logscale x
set logscale y
set xrange [*:1]
set yrange [*:1]

#set dgrid3d 1000,1000
set hidden3d
set isosamples 100,100

set ylabel "y"
set xlabel "x"
set zlabel "corr"
tau=0.1
cutoff=1e-5;
f(x,y)=(x>cutoff && x<1 && y>cutoff&&y<1) ? 1/(2*pi)* sqrt(x*y)* 1/sqrt(x*y*(1-x-y)) * (exp((-pi*tau*tau)/(1-x-y)) - exp((-4*pi*tau*tau)/(1-x-y))):0/0; #Nothing if sum above 1

splot "./Results/histogram2d_simple.dat" u 1:2:3 t"numerical" pt 0 with points,\
f(x,y) t"theory" lt 1 lc rgb"red" with lines




