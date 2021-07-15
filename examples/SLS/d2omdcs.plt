reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set xrange [0:30]
set style line 1 lt 1 lc 1
plot \
'SLS39631.mdcs01' u 1:((column(2)+0.0)+0.0) notitle w l ls 1
