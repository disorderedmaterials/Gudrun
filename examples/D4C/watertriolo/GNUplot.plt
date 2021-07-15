reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set logscale y
set style line 1 lt 1ps 0 lc 1
set style line 2 lt 1ps 0 lc 4
set style line 3 lt 1ps 0 lc 7
x=0
y=0
i=-1
plot \
'178834.mint01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+1.0)*1.0):1/0):1/0)) title 'D2O D4C' w l ls 1, \
'178745.mint01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+1.0)*2.0):1/0):1/0)) title 'H2O D4C' w l ls 2, \
'178724.mint01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+1.0)*3.0):1/0):1/0)) title 'Ni Powder' w l ls 3
