set terminal epslatex size 15cm,7.5cm font ',10' color
#input header '\small'


set key samplen 3 box height .5 top right
set ylabel 'Signal /[A.u.]'
set xlabel '$t$ /[$\unit{ms}$]'
set grid
set format y '%g'


A2 =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $1}'")+0
p  =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $2}'")+0
R1 =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $3}'")+0
R2 =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $4}'")+0

A1 =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $5}'")+0
R  =system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $6}'")+0

bkg=system("head -1 fit-param_1_0ms_xwpp.tsv | awk '{print $7}'")+0


s='\footnotesize '

set grid

set key width -3.6
set output "waveform.tex"
plot [-10:200] [-0.1:1.5] "1_0ms_xwpp_WAVEFORM.tsv" using 1:($2-bkg) with lines lt -1 lc 3 title s.'Wave-form',\
                     "1_0ms_xwpp_PEAKS.tsv" using 1:($2-bkg) lt 1 pt 6 ps 2 lc 7 title s.'Peaks',\
                     A2*( p*exp(-R1*x) + (1-p)*exp(-R2*x) )  dt 1 lw 4 lc 0 title s.'Dual mode fit',\
                     A1*exp(-R*x)                            dt 2 lw 4 lc 0 title s.'Single mode fit'



#set terminal epslatex size 8cm,7.5cm font ',10' color
set format y '$10^%L$'
set logscale y
set output "peaks.tex"
plot [-10:200] [5e-2:2] "1_0ms_xwpp_PEAKS.tsv" using 1:($2-bkg) pt 6 ps 1.2 lc 7 title s.'Peaks',\
                        A2*( p*exp(-R1*x) + (1-p)*exp(-R2*x) ) dt 1 lw 4 lc 0 title s.'Dual mode fit',\
                        A1*exp(-R*x)                           dt 2 lw 4 lc -1 title s.'Single mode fit'


set out
