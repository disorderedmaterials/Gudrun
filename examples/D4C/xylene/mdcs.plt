reset
set xlabel 'Q [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set title 'mint 20C'
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
x=0
y=0
i=-1
plot \
'225923.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+1.0)+0.0):1/0):1/0)) notitle w l ls 1, \
'225923.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(4),(i>0?(lasty!=0?((lasty+1.0)+0.0):1/0):1/0)) notitle w l ls 2, \
'225923.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(6),(i>0?(lasty!=0?((lasty+1.0)+0.0):1/0):1/0)) notitle w l ls 3, \
'225923.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(8),(i>0?(lasty!=0?((lasty+1.0)+0.0):1/0):1/0)) notitle w l ls 4, \
'225923.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(10),(i>0?(lasty!=0?((lasty+1.0)+0.0):1/0):1/0)) notitle w l ls 5, \
'225923.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(12),(i>0?(lasty!=0?((lasty+1.0)+0.0):1/0):1/0)) notitle w l ls 6, \
'225923.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(14),(i>0?(lasty!=0?((lasty+1.0)+0.0):1/0):1/0)) notitle w l ls 7, \
'225923.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(16),(i>0?(lasty!=0?((lasty+1.0)+0.0):1/0):1/0)) notitle w l ls 8, \
'225923.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(18),(i>0?(lasty!=0?((lasty+1.0)+0.0):1/0):1/0)) notitle w l ls 9, \
'225923.dcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(20),(i>0?(lasty!=0?((lasty+1.0)+0.0):1/0):1/0)) notitle w l ls 10
