#script mc
#set term png
#set output 'fig5E_mc.png'

set title 'Energy Error' font ',18'

set xtics font ',11'
set ytics font ',11'

set xrange [0:1000000]

set xlabel 'm' font ',18'
set ylabel '{/Symbol s}' font ',18'

plot "R_binn2.00.txt" u ($1):($3) t 'T=2.0' w l,\
