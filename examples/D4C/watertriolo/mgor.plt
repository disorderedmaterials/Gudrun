reset
set terminal postscript eps enhanced colour "Times-Roman" 22 
set encoding iso_8859_1
set output 'mint.eps'
set xlabel 'r [1/{\305}]'
set ylabel 'Diff. CS [barns/sr/atom]'
set xrange [0:10]
set style line 1 lt 1ps 0 lc 1
set style line 2 lt 1ps 0 lc 3
set style line 3 lt 1ps 0 lc 4
set style line 4 lt 1ps 0 lc 5
x=0
y=0
i=-1
plot \
'178834.mgor01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.5)+0.0):1/0):1/0)) title 'D2O D4C' w l ls 1, \
'/home/aks45/Gudrun/Gudrun2014/Gudrun/run/SLS/SLS39631.mgor01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.5)+0.0):1/0):1/0)) title 'D2O SANDALS' w l ls 2, \
'178745.mgor01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.5)+-0.5):1/0):1/0)) title 'H2O D4C' w l ls 3, \
'/home/aks45/Gudrun/Gudrun2014/Gudrun/run/SLS/SLS39633.mgor01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.5)+-0.5):1/0):1/0)) title 'H2O SANDALS' w l ls 4
set output
