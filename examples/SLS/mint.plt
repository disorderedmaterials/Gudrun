reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set xrange [0:30]
set style line 1 lt 1ps 0 lc 1
set style line 2 lt 1ps 0 lc 2
set style line 3 lt 1ps 0 lc 3
plot \
'SLS39631.mint01' u 1:((column(2)+0.0)+0.0) notitle w l ls 1, \
'SLS39632.mint01' u 1:((column(2)+0.0)+0.0) notitle w l ls 2, \
'SLS39633.mint01' u 1:((column(2)+0.0)+0.0) notitle w l ls 3
