reset
set xlabel 'Q [1/A]'
set ylabel 'Diff. CS'
set style line 1 lt 1 lc 1
plot \
'090520d_H2O.int01' u 1:((column(2)+0.0)+0.0) notitle w l ls 1
