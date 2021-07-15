reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set xrange [0:30]
set style line 1 lt 1ps 0 lc 1
set style line 2 lt 1ps 0 lc 2
set style line 3 lt 1ps 0 lc 3
x=0
y=0
i=-1
plot \
'SLS50875.mdcs01' u 1:((column(2)+0.0)+0.0) notitle w l ls 1, \
'SLS50876.mdcs01' u 1:((column(2)+0.0)+0.0) notitle w l ls 2, \
'SLS50994.mdcs01' u 1:((column(2)+0.0)+0.0) notitle w l ls 3
