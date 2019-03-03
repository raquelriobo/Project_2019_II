#script binning
#set term postscript
#set output 'Epotbin.eps'
set term png
set output 'Epotbin.png'

set title 'Potential Energy Error' font ',18'

set xtics font ',11'
set ytics font ',11'

set xrange [0:2200]
set yrange [0:0.52]

set xlabel 'm' font ',18'
set ylabel 's' font ',18'

set key bottom inside

plot "Epotbinning.txt" u 1:2 t 'Potential Energy' w lp pt 3 lc rgb "purple"

#set term postscript
#set output 'Ekinbin.eps'
set term png
set output 'kinbin.png'
set title 'Kinetic Energy Error' font ',18'

set xtics font ',11'
set ytics font ',11'

set xrange [0:2200]
set yrange [0:0.52]

set xlabel 'm' font ',18'
set ylabel 's' font ',18'

set key bottom inside

plot "Ekinbinning.txt" u 1:2 t 'Kinetic Energy' w lp pt 3 lc rgb "blue"

#set term postscript
#set output 'Etotbin.eps'
set term png
set output 'totbin.png'
set title 'Total Energy Error' font ',18'

set xtics font ',11'
set ytics font ',11'

set xrange [0:2200]
set yrange [0:0.055555]

set xlabel 'm' font ',18'
set ylabel 's' font ',18'

set key bottom inside

plot "Etotbinning.txt" u 1:2 t 'Total Energy' w lp pt 3 lc rgb "red"

#set term postscript
#set output 'Etotbin.eps'
set term png
set output 'comparebin.png'
set title 'Compare Energy Error' font ',18'

set xtics font ',11'
set ytics font ',11'

set xrange [0:2200]
set yrange [0:0.52]

set xlabel 'm' font ',18'
set ylabel 's' font ',18'

set key center right inside

plot "Etotbinning.txt" u 1:2 t 'Total Energy' w lp pt 3 lc rgb "red",\
"Ekinbinning.txt" u 1:2 t 'Kinetic Energy' w lp pt 3 lc rgb "blue",\
"Epotbinning.txt" u 1:2 t 'Potential Energy' w lp pt 3 lc rgb "purple"
