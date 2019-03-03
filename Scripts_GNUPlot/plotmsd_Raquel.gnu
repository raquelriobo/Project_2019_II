#script dm
set term png
set output 'MSD_PF2.png'
set xtics font ',11'
set ytics font ',11'
set xlabel 'time (ps)' font ',18'
set ylabel 'msd (Angstrom^{2})' font ',18'
set key box nobox bmargin inside bottom
f(x)=a+b*x
fit f(x) "v_verlet_2_0.2_msd.txt" u ($1):($2) via a,b
plot "v_verlet_2_0.2_msd.txt" u ($1):($2) t 'Mean Square Displacement' w l,\
      f(x) ls 4 t 'f(x)=mx+n'
