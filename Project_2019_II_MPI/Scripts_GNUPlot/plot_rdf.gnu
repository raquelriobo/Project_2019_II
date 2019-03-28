#Script plot RDF
set term pngcairo

set term postscript enhanced color
set output 'rdf.eps'
set xlabel 'r/sigma'
set ylabel 'g(r)'
set title 'Radial distribution function'
unset key

plot 'radial.txt' u 1:2 w l
