#set xrange [0:10]
#set yrange [0:300]
set ylabel "Tiempo por iteración (s)"
set xlabel "Número de procesadores"
unset key
set tics nomirror
set terminal epslatex size 3.2,2.4
set output "timef_iter_3.tex"
plot "total_3.dat" u 1:5 w l lw 2 lc rgb "#1e4ab8"
