reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set xrange [0.5:60]
set logscale x
set style line 1 lt 1 lc 2
plot \
'GEM46509.mint01' u 1:(column(2)+0.0) notitle w l ls 1
