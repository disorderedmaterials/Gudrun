reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set style line 1 lt 1ps 0 lc 1
set style line 2 lt 1ps 0 lc 2
set style line 3 lt 1ps 0 lc 3
set style line 4 lt 1ps 0 lc 4
set style line 5 lt 1ps 0 lc 5
set style line 6 lt 1ps 0 lc 6
set style line 7 lt 1ps 0 lc 7
set style line 8 lt 1ps 0 lc 8
set style line 9 lt 1ps 0 lc 9
set style line 10 lt 1ps 0 lc 10
set style line 11 lt 1ps 0 lc 11
set style line 12 lt 1ps 0 lc 12
set style line 13 lt 1ps 0 lc 13
set style line 14 lt 1ps 0 lc 14
set style line 15 lt 1ps 0 lc 15
set style line 16 lt 1ps 0 lc 16
set style line 17 lt 1ps 0 lc 17
set style line 18 lt 1ps 0 lc 18
plot \
'SLS39631.sub01' u 1:((column(2)+0.0)+0.0) notitle w l ls 1, \
'SLS39631.sub01' u 1:((column(4)+0.0)+0.0) notitle w l ls 2, \
'SLS39631.sub01' u 1:((column(6)+0.0)+0.0) notitle w l ls 3, \
'SLS39631.sub01' u 1:((column(8)+0.0)+0.0) notitle w l ls 4, \
'SLS39631.sub01' u 1:((column(10)+0.0)+0.0) notitle w l ls 5, \
'SLS39631.sub01' u 1:((column(12)+0.0)+0.0) notitle w l ls 6, \
'SLS39631.sub01' u 1:((column(14)+0.0)+0.0) notitle w l ls 7, \
'SLS39631.sub01' u 1:((column(16)+0.0)+0.0) notitle w l ls 8, \
'SLS39631.sub01' u 1:((column(18)+0.0)+0.0) notitle w l ls 9, \
'SLS39631.sub01' u 1:((column(20)+0.0)+0.0) notitle w l ls 10, \
'SLS39631.sub01' u 1:((column(22)+0.0)+0.0) notitle w l ls 11, \
'SLS39631.sub01' u 1:((column(24)+0.0)+0.0) notitle w l ls 12, \
'SLS39631.sub01' u 1:((column(26)+0.0)+0.0) notitle w l ls 13, \
'SLS39631.sub01' u 1:((column(28)+0.0)+0.0) notitle w l ls 14, \
'SLS39631.sub01' u 1:((column(30)+0.0)+0.0) notitle w l ls 15, \
'SLS39631.sub01' u 1:((column(32)+0.0)+0.0) notitle w l ls 16, \
'SLS39631.sub01' u 1:((column(34)+0.0)+0.0) notitle w l ls 17, \
'SLS39631.sub01' u 1:((column(36)+0.0)+0.0) notitle w l ls 18
