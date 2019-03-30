#GNUPlot script to generate energies, pressure and temperature figures

set term postscript enhanced color 
set output 'energies.eps'
set xlabel 'Time (ps)' 
set ylabel 'Energy (KJ/mol)'
set title 'Energy vs time'
set key outside

plot 'E_kinpot.txt' u 1:2 t 'Potential Energy' w l lc rgb 'web-green',\
     '' u 1:3 t 'Kinetic Energy' w l lc rgb 'web-blue',\
     'E_tot.txt' u 1:2 t 'Total Energy' w l lc rgb 'orange'

set output 'pressure.eps'
set ylabel 'Pressure (Pa)'
set title 'Pressure vs time'
plot 'PressTemp.txt'  u 1:2  w l notitle lc rgb 'pink'

set output 'temperature.eps'
set ylabel 'Temperature (K)'
set title 'Temperature vs time'
plot 'PressTemp.txt'  u 1:3 w l notitle  lc rgb 'red'

