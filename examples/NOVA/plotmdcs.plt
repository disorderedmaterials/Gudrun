reset
set xlabel 'Q [1/A]'
set ylabel 'Diff CS [barn/sr/atom]'
set title 'mdcs plot'
set xrange [0.01:50]
set logscale x
set logscale y
set style line 1 lt 1ps 0 lc 1
set style line 2 lt 1ps 0 lc 2
set style line 3 lt 1ps 0 lc 3
set style line 4 lt 1ps 0 lc 4
set style line 5 lt 1ps 0 lc 5
set style line 6 lt 1ps 0 lc 7
set style line 7 lt 1ps 0 lc 8
x=0
y=0
i=-1
plot \
'014744_10steps.mdcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.0)*1.0):1/0):1/0)) title 'PSD8' w l ls 1, \
'014775_10steps.mdcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.0)*1.0):1/0):1/0)) title 'PSD5' w l ls 2, \
'014754_10steps.mdcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.0)*1.0):1/0):1/0)) title 'PSH8' w l ls 3, \
'014764_10steps.mdcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.0)*1.0):1/0):1/0)) title 'PSD8D5' w l ls 4, \
'014733_10steps.mdcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.0)*1.0):1/0):1/0)) title 'PSH8D8' w l ls 5, \
'014885_10steps.mdcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.0)*1.0):1/0):1/0)) title 'PSH8D5' w l ls 6, \
'014765_10steps.mdcs01' u (xval=(column(1)+x)*0.5,x=column(1),xval):(i=i+1,lasty=y,y=column(2),(i>0?(lasty!=0?((lasty+0.0)*1.0):1/0):1/0)) title 'Al Container' w l ls 7
