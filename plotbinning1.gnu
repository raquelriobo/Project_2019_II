#script mc
#set term png
#set output 'BinningL20.eps'

set xlabel 'm' font ',18'
set ylabel '{/Symbol s}' font ',18'

set xtics font ',11'
set ytics font ',11'

set key box nobox bottom inside

a1=0.002
a2=0.002

b1=0.005
b2=0.005

tau1=1
tau2=1

fit1(x)=a1-b1*exp(-x/tau1)
fit2(x)=a2-b2*exp(-x/tau2)
fit fit1(x) 'Temp2.0L20.txt' u 1:2 via a1,b1,tau1
fit fit2(x) 'Temp4.0L20.txt' u 1:2 via a2,b2,tau2
plot 'Temp2.0L20.txt' w p pt 4 lc rgb 'red' t 'T=2',\
     'Temp4.0L20.txt' w p pt 4 ps 1 lc rgb 'web-blue' t 'T=4',\
     fit1(x) w l lw 1 lt 1 lc rgb 'red' t 'Fit T=2',\
     fit2(x) w l lw 1 lt 1 lc rgb 'blue' t 'Fit T=4'
