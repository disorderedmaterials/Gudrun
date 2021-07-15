reset
set xlabel 'r [A]'
set ylabel 'g(r)'
set title '[barns/sr/atom]'
set xrange [0:12]
set style line 1 lt 1ps 0 lc 1
plot \
'SLS39633.mgor01' u 1:((column(2)+0.0)+0.0) notitle w l ls 1
