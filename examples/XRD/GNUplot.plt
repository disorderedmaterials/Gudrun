reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set style line 1 lt 1ps 0 lc 1
set style line 2 pt 1 lc 2
x=0
y=0
i=-1
plot \
'090523a_SiO2.int01' u 1:((column(2)+0.0)+0.0) notitle w l ls 1, \
'090520d_H2O.int01' u 1:((column(2)+0.0)+0.0) notitle w p ls 2
