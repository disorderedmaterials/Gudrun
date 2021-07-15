reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set style line 1 lt 1ps 0 lc 1
set style line 2 lt 1ps 0 lc 2
set style line 3 lt 1ps 0 lc 3
set style line 4 lt 1ps 0 lc 4
set style line 5 pt 1 lc 5
x=0
y=0
i=-1
plot \
'008266.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.0)+0.0):1/0):1/0)) notitle w l ls 1, \
'008266.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(4),(i>0?(lasty!=0?((lasty+0.0)+0.0):1/0):1/0)) notitle w l ls 2, \
'008266.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(6),(i>0?(lasty!=0?((lasty+0.0)+0.0):1/0):1/0)) notitle w l ls 3, \
'008266.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(8),(i>0?(lasty!=0?((lasty+0.0)+0.0):1/0):1/0)) notitle w l ls 4, \
'008266.mdcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.0)+0.0):1/0):1/0)) notitle w p ls 5
