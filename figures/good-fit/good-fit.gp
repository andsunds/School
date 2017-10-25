set terminal epslatex size 15cm,7.5cm font ',10' color
#input header '\small'


set key samplen 2 box height .5 top left
set ylabel '$1-R^2$'
set xlabel '$\tau$ /[$\unit{ms}$]'
set grid
set format y '$10^%L$'
set logscale xy

s='\footnotesize '

set grid

set key width -3.6
set output "fit_Rsq.tex"
plot [3e-1:30] [1e-4:1e0] "d-water_Rsq.tsv" using 1:(1-$2) lt 1 pt 9 ps 1.2 lc rgb 'blue' title s.'Liquid water',\
                          "d-water_Rsq.tsv" using 1:(1-$3) lt 1 pt 8 ps 1.2
lc rgb 'blue' title '',\
                       "wet-pp_Rsq.tsv" using 1:(1-$2) lt 1 pt 5 ps 1.2 lc rgb 'green' title s.'Wet paper',\
                       "wet-pp_Rsq.tsv" using 1:(1-$3) lt 1 pt 4 ps 1.2 lc rgb 'green' title '',\
                       "ex-wet-pp_Rsq.tsv" using 1:(1-$2) lt 1 pt 7 ps 1.2 lc rgb 'red' title s.'Ex. wet paper',\
                       "ex-wet-pp_Rsq.tsv" using 1:(1-$3) lt 1 pt 6 ps 1.2 lc rgb 'red' title '',\



set out
