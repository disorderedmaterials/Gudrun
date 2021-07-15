reset
set xlabel 'r [A]'
set ylabel 'g(r)'
set style line 1 lt 1 lc 1
plot \
'090520d_H2O.gofr' u 1:((column(2)+0.0)+0.0) notitle w l ls 1
