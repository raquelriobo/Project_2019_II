#set xrange [0:10]
set yrange [0:1.5]
set ylabel "Eficiencia"
set xlabel "Número de procesadores"
#unset key
#set key bottom
set tics nomirror
#set logscale
set terminal epslatex size 3.2,2.4
set output "eficiencia_total_3.tex"
plot "total_3.dat" u 1:((71.2950473/($2))/($1)) w l lw 2 lc rgb "#1e4ab8" t "Eficiencia", 1 w l lw 2 lc rgb "red" t "Ideal"

