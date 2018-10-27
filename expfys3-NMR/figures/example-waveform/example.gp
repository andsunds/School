set terminal epslatex size 15cm,9cm font ',10' color
#input header '\small'


set key box height .5 top right
set ylabel 'Oscilloscope signal /[V]'
set xlabel '$t$ /[$\unit{ms}$]'
set grid
set format y '%g'

set grid

s='\footnotesize '

A28  =system("head -1 fit-param_28ms_dw.tsv | awk '{print $1}'")+0
R28  =system("head -1 fit-param_28ms_dw.tsv | awk '{print $2}'")+0
bkg28=system("head -1 fit-param_28ms_dw.tsv | awk '{print $3}'")+0

A50  =system("head -1 fit-param_50ms_dw.tsv | awk '{print $1}'")+0
R50  =system("head -1 fit-param_50ms_dw.tsv | awk '{print $2}'")+0
bkg50=system("head -1 fit-param_50ms_dw.tsv | awk '{print $3}'")+0

set key width -2.2
set output "waveform.tex"
plot [-10:1000] [-0.1:2] "28ms_dw_WAVEFORM.tsv" using 1:2 with lines lt -1 lc 3 title s.'Waveform, $\tau=\unit[28]{ms}$',\
                         "50ms_dw_WAVEFORM.tsv" using 1:2 with lines lt -1 lc 2 title s.'Waveform, $\tau=\unit[50]{ms}$',\
                         "28ms_dw_PEAKS.tsv" using 1:2 pt 12 ps 2 lc 7 title s.'Peaks, $\tau=\unit[28]{ms}$',\
                         "50ms_dw_PEAKS.tsv" using 1:2 pt 6  ps 2 lc 7 title s.'Peaks, $\tau=\unit[50]{ms}$',\
                         A28*exp(-R28*x)+bkg28 dt 4 lw 4 lc 0 title s.'Fit, $\tau=\unit[28]{ms}$',\
                         A50*exp(-R50*x)+bkg50 dt 1 lw 4 lc 0 title s.'Fit, $\tau=\unit[50]{ms}$'




#set terminal epslatex size 8cm,7.5cm font ',10' color
set key bottom right#at 380,2.5e-1#
set ylabel 'Peaks /[V] (above noise background level)'
set format y '$10^%L$'
set logscale y

A2 =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $1}'")+0
p  =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $2}'")+0
R1 =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $3}'")+0
R2 =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $4}'")+0

A1 =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $5}'")+0
R  =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $6}'")+0

bkg=system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $7}'")+0

set key width -7
set output "peaks.tex"
plot [-10:450] [7e-2:2] "1_0ms_xwpp_PEAKS.tsv" using 1:($2-bkg) pt 6 ps 1.2 lc 7 title s.'Peaks, ex. wet paper',\
                        A2*( p*exp(-R1*x) + (1-p)*exp(-R2*x) ) dt 1 lw 4 lc 0 title s.'Dual mode fit',\
                        A1*exp(-R*x)                           dt 2 lw 4 lc -1 title s.'Single mode fit',\
                        "28ms_dw_PEAKS.tsv" using 1:($2-bkg50) pt 12 ps 2 lc 7 title s.'Peaks, liquid water',\
                        A28*exp(-R28*x) dt 4 lw 4 lc 0 title s.'Single mode fit '


set out
