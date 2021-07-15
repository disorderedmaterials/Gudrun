reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set style line 1 lt 1ps 0 lc 1
set style line 2 lt 1ps 0 lc 2
set style line 3 lt 1ps 0 lc 3
set style line 4 lt 1ps 0 lc 4
set style line 5 pt 1 lc 5
plot \
'008162.dcs01' u 1:((column(2)+0.0)+0.0) notitle w l ls 1, \
'008162.dcs01' u 1:((column(4)+0.0)+0.0) notitle w l ls 2, \
'008162.dcs01' u 1:((column(6)+0.0)+0.0) notitle w l ls 3, \
'008162.dcs01' u 1:((column(8)+0.0)+0.0) notitle w l ls 4, \
'008162.mdcs01' u 1:((column(2)+0.0)+0.0) notitle w p ls 5
