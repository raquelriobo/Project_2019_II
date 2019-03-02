#script dm
set term png
set output 'momentum.png'
set xtics font ',11'
set ytics font ',11'
set xlabel 'Time' font ',18'
set ylabel 'Momentum' font ',18'
set xrange [0:3]
set key box

plot 'v_verlet_momnt.txt' u 1:2 t 'x' w l lc rgb 'pink',\
      '' u 1:3 t 'y' w l lc rgb 'skyblue',\
      '' u 1:4 t 'z' w l lc rgb 'orange'
