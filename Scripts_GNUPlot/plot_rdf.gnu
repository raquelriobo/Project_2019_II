#Script plot RDF
set term pngcairo
set output 'rdf.png'
set xtics font ',11'
set ytics font ',11'
set xlabel 'r/sigma' font ',18'
set ylabel 'g(r)' font ',18'
#set xrange [0:3]
set key box outside nobox

plot 'radial.txt' u 1:2 w lp
