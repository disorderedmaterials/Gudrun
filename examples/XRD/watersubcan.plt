reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set style line 1 lt 1 lc 1
set style line 2 lt 1 lc 2
plot \
'090520d_H2O.subcan' u 1:((column(2)+0.0)+0.0) notitle w l ls 1, \
'090520d_H2O.subcan' u 1:((column(4)+0.0)+0.0) notitle w l ls 2
