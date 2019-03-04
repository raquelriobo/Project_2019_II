#script gnuplot binning
set term postscript color enhanced
set output 'Epotbin.eps'

set title 'Potential Energy Error'
set xlabel 'm'
set ylabel '{/Symbol s}'

unset key

plot "Epotbinning.txt" u 1:2 t 'Potential Energy' w lp pt 3 lc rgb "purple"

set output 'Ekinbin.eps'
set title 'Kinetic Energy Error'

plot "Ekinbinning.txt" u 1:2 t 'Kinetic Energy' w lp pt 3 lc rgb "blue"

set output 'Etotbin.eps'
set title 'Total Energy Error' 

plot "Etotbinning.txt" u 1:2 t 'Total Energy' w lp pt 3 lc rgb "red"

set output 'Tempbin.eps'
set title 'Temperature Error'

plot "Tempbinning.txt" u 1:2 w lp pt 3 lc rgb "blue"

set output 'Pressurebin.eps'
set title 'Pressure Error'

plot "Presbinning.txt" u 1:2 w lp pt 3 lc rgb "blue"
