reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set xrange [0.:30]
set style line 1 lt 1ps 0 lc 1
x=0
y=0
i=-1
plot \
'GEM46509.mint01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.0)+0.0):1/0):1/0)) title 'GEM46509.mint01' w l ls 1
