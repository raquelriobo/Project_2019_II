set term postscript enhanced color
set output 'total_momentum.eps'
set xlabel 'Time (reduced units)'
set ylabel 'Momentum (reduced units)'
set title 'Total momentum vs time'
set box outside
plot 'Momentum.txt' u 1:2 t 'Px' w l lc rgb 'pink',\
      '' u 1:3 t 'Py' w l lc rgb 'skyblue',\
      '' u 1:4 t 'Pz' w l lc rgb 'orange'
