#script dm
set term png
set output 'results.png'
set xtics font ',11'
set ytics font ',11'
set xlabel 'Time' font ',18'
set ylabel 'Energy' font ',18'
set xrange [0:3]
set key box outside nobox

plot 'Results.txt' u 1:2 t 'potential' w l lc rgb 'web-green',\
      '' u 1:3 t 'kinetic' w l lc rgb 'web-blue',\
      '' u 1:4 t 'total' w l lc rgb 'orange'
      '' u 1:5 t 'pressure' w l lc rgb 'pink'
      '' u 1:6 t 'temperature' w l lc rgb 'red'
