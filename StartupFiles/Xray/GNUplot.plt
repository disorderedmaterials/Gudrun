reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set style line 1 lt 1ps 0 lc 1
x=0
y=0
i=-1
plot \
'WO3-corrected.XRDtxt' u 1:((column(2)+0.0)+0.0) notitle w l ls 1
