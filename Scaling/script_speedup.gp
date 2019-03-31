#set xrange [0:10]
set yrange [-1:100]
set ylabel "Speedup"
set xlabel "Número de procesadores"
#unset key
set key bottom
set tics nomirror
set logscale
set terminal epslatex size 3.2,2.4
set output "speedup_total_3.tex"
plot "total_3.dat" u 1:(71.2950473/($2)) w l lw 2 lc rgb "#1e4ab8" t "Speedup", x w l lw 2 lc rgb "red" t "Ideal"

