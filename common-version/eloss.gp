# GS version of the energy loss plot with quite a larger statistics

reset
set term pdfcairo enhanced
set out 'eloss.pdf'
set colors classic

fn='res/eloss-tmax_3.0-xmin_1e-4-eps_1e-9.res'
fn2='res/eloss-tmax_3.0-xmin_1e-3-eps_1e-9.res'
fn3='res/eloss-tmax_3.0-xmin_1e-3-eps_1e-8.res'

fndef=fn

set label 1 '{/Symbol e}=10^{-9}, x_{min}=10^{-4}' at graph 0.05,0.07
xmin=0.0001

set xlabel '{/Symbol t}'

set ylabel '1-(energy loss)'
set log y
set yrange [1e-8:1]
set format y '10^{%T}'
set ytics add ("1" 1 0)
set ytics add ("0.1" 0.1 0)



nev=system('grep "nev = " '.fndef.' | awk "{print \$5}"')+0

g(t)=exp(-a*t*t-b*t)*(1.0-2.0*t*sqrt(xmin))
ga(t,a,b,xm)=exp(-a*t*t-b*t)*(1.0-2.0*t*sqrt(xm))

plot fndef u 1:4 w linesp dt 1 lc 3 lw 2 pt 7 ps 0.6 t 'simple (numeric)',\
     fndef u 1:($4-sqrt($5/(1.0*nev))) w l dt 1 lc 3 lw 1 not,\
     fndef u 1:($4+sqrt($5/(1.0*nev))) w l dt 1 lc 3 lw 1 not,\
     ga(x,pi,0.0,xmin)                                  w l      dt 2 lc 3 lw 3             t 'simple (analytic)',\
     fndef u 1:6 w linesp dt 1 lc 1 lw 2 pt 9 ps 0.6 t 'full (numeric)'

# try to extract the propoerties by fitting a Gaussian
set fit errors
fit [0.0:1.85] log(g(x)) fndef u 1:(log($4)) via a,b
a_simple=a
a_simple_err=a_err
b_simple=b
b_simple_err=b_err
fit [0.0:2.25] log(g(x)) fndef u 1:(log($6)) via a,b
a_full=a
a_full_err=a_err
b_full=b
b_full_err=b_err

plot fndef u 1:4 w linesp dt 1 lc 3 lw 2 pt 7 ps 0.6 t 'simple (numeric)',\
     ga(x,pi,0.0,xmin)                                 w l      dt 2 lc 3 lw 3             t 'simple (analytic)',\
     fndef u 1:6 w linesp dt 1 lc 1 lw 1 pt 9 ps 0.6 t 'full (numeric)',\
     ga(x,a_simple,b_simple,xmin)                      w l      dt 1 lc 7 lw 2             t 'fit',\
     ga(x,a_full,b_full,xmin)                          w l      dt 1 lc 7 lw 1 not


print "Simple : a=",a_simple," +- ",a_simple_err," (expected pi, i.e. pull=",(a_simple-pi)/a_simple_err,")"
print "         b=",b_simple," +- ",b_simple_err," (expected  0, i.e. pull=",b_simple/b_simple_err,")"
print "Full   : a=",a_full," +- ",a_full_err
print "         b=",b_full," +- ",b_full_err

set ylabel 'variance'
plot fndef u 1:5 w linesp dt 1 lc 3 lw 2 pt 7 ps 0.6 t 'simple (numeric)',\
     fndef u 1:3 w l      dt 2 lc 3 lw 3             t 'simple (analytic)',\
     fndef u 1:7 w linesp dt 1 lc 1 lw 2 pt 9 ps 0.6 t 'full (numeric)'

set yrange [1e-12:1]
plot fndef u 1:4 w l dt 1 lc 3 lw 2 not,\
     for [i=1:30] 'res/parts/eloss-tmax_3.0-xmin_1e-4-eps_1e-9-rseq'.sprintf("%d",i).'.res' u 1:4 w l dt 1 lc 7 lw 0.5 not

unset log y
set ylabel 'ratio to th expected'
set format y "%g"
set yrange [0:1.2]
plot fndef u 1:(($4-sqrt($5/(1.0*nev)))/ga($1,pi,0.0,xmin)):(($4+sqrt($5/(1.0*nev)))/ga($1,pi,0.0,xmin)) w filledcurve fs solid 0.3 lc 3 not,\
     fndef u 1:($4/ga($1,pi,0.0,xmin))                 w linesp dt 1 lc 3 lw 2 pt 7 ps 0.6 t 'numeric',\
     ga(x,a_simple,b_simple,xmin)/ga(x,pi,0.0,xmin)    w l      dt 1 lc 7 lw 2             t 'fit'
     

set ylabel 'ratio to fit'
plot fndef u 1:(($6-sqrt($7/(1.0*nev)))/ga($1,a_full,b_full,xmin)):(($6+sqrt($7/(1.0*nev)))/ga($1,a_full,b_full,xmin)) w filledcurve fs solid 0.3 lc 3 not,\
     fndef u 1:($6/ga($1,a_full,b_full,xmin))                      w linesp dt 1 lc 3 lw 2 pt 7 ps 0.6 t 'numeric/fit'

unset label 1
set ylabel 'ratio to th expectation'
set key bottom left reverse Left
plot 1.0 w l dt 2 lc 7 lw 1 not,\
     fn3 u 1:(($4-sqrt($5/(1.0*nev)))/ga($1,pi,0.0,0.001)):(($4+sqrt($5/(1.0*nev)))/ga($1,pi,0.0,0.001)) w filledcurve fs solid 0.3 transparent lc 3 not,\
     fn3 u 1:($4/ga($1,pi,0.0,0.001))                      w linesp dt 1 lc 3 lw 2 pt 7 ps 0.6 t '{/Symbol e}=10^{-8}, x_{min}=10^{-3}',\
     fn2 u 1:(($4-sqrt($5/(1.0*nev)))/ga($1,pi,0.0,0.001)):(($4+sqrt($5/(1.0*nev)))/ga($1,pi,0.0,0.001)) w filledcurve fs solid 0.3 transparent lc 1 not,\
     fn2 u 1:($4/ga($1,pi,0.0,0.001))                      w linesp dt 1 lc 1 lw 2 pt 6 ps 0.6 t '{/Symbol e}=10^{-9}, x_{min}=10^{-3}',\
     fn  u 1:(($4-sqrt($5/(1.0*nev)))/ga($1,pi,0.0,0.0001)):(($4+sqrt($5/(1.0*nev)))/ga($1,pi,0.0,0.0001)) w filledcurve fs solid 0.2 transparent lc 7 not,\
     fn  u 1:($4/ga($1,pi,0.0,0.0001))                     w linesp dt 1 lc 7 lw 2 pt 8 ps 0.6 t '{/Symbol e}=10^{-9}, x_{min}=10^{-4}'

set ylabel 'ratio of E_{kept} full/simple'
set key bottom top right noreverse Right
set yrange [0:10]
set grid
plot 1.0 w l dt 2 lc 7 lw 1 not,\
     fn3 u 1:($6/$4)  w linesp dt 1 lc 3 lw 2 pt 7 ps 0.6 t '{/Symbol e}=10^{-8}, x_{min}=10^{-3}',\
     fn2 u 1:($6/$4)  w linesp dt 1 lc 1 lw 2 pt 6 ps 0.6 t '{/Symbol e}=10^{-9}, x_{min}=10^{-3}',\
     fn  u 1:($6/$4)  w linesp dt 1 lc 7 lw 2 pt 8 ps 0.6 t '{/Symbol e}=10^{-9}, x_{min}=10^{-4}'

set ylabel 'ratio of E_{loss} full/simple'
set key bottom top right noreverse Right
set yrange [0:1.2]
set grid
plot 1.0 w l dt 2 lc 7 lw 1 not,\
     fn3 u 1:((1-$6)/(1-$4))  w linesp dt 1 lc 3 lw 2 pt 7 ps 0.6 t '{/Symbol e}=10^{-8}, x_{min}=10^{-3}',\
     fn2 u 1:((1-$6)/(1-$4))  w linesp dt 1 lc 1 lw 2 pt 6 ps 0.6 t '{/Symbol e}=10^{-9}, x_{min}=10^{-3}',\
     fn  u 1:((1-$6)/(1-$4))  w linesp dt 1 lc 7 lw 2 pt 8 ps 0.6 t '{/Symbol e}=10^{-9}, x_{min}=10^{-4}'

set ylabel 'ratio of fluctuations to loss'
set key bottom top right noreverse Right
set yrange [0:0.6]
set grid
plot 1.0/sqrt(3.0) w l dt 2 lc 7 lw 1 not,\
     fn u 1:(sqrt($5)/(1-$4))  w linesp dt 1 lc 1 lw 2 pt 6 ps 0.6 t 'simple',\
     fn u 1:(sqrt($7)/(1-$6))  w linesp dt 1 lc 7 lw 2 pt 8 ps 0.6 t 'full'


set out
