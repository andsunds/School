set terminal epslatex size 15cm,9cm font ',10' color
#input header '\small'


set key samplen 2 box height .5 top left
set ylabel '$1-r^2$'
set xlabel '$\tau^2$ /[$\unit{ms^2}$]'
set grid
set format y '$10^%L$'
set logscale xy

s='\footnotesize '

set grid

set key width -6
set output "fit_Rsq.tex"
plot [1e-1:1e3] [1e-4:1e0] "d-water_Rsq.tsv" using ($1**2):(1-$2) lt 1 pt 9 ps 1.2 lc 2 title s.'Liquid water',\
                          "d-water_Rsq.tsv" using ($1**2):(1-$3) lt 1 pt 8 ps 1.2 lc 2 title '',\
                          "ex-wet-pp_Rsq.tsv" using ($1**2):(1-$2) lt 1 pt 7 ps 1.5 lc 7 title s.'Ex. wet paper',\
                          "ex-wet-pp_Rsq.tsv" using ($1**2):(1-$3) lt 1 pt 6 ps 1.5 lc 7 title '',\
                          "wet-tp_Rsq.tsv" using ($1**2):(1-$2) lt 1 pt 13 ps 1.5 lc 6 title s.'Wet tissue paper',\
                          "wet-tp_Rsq.tsv" using ($1**2):(1-$3) lt 1 pt 12 ps 1.5 lc 6 title '',\
                          #"wet-pp_Rsq.tsv" using 1:(1-$2) lt 1 pt 5 ps 1.2 lc 2 title s.'Wet paper',\
                          #"wet-pp_Rsq.tsv" using 1:(1-$3) lt 1 pt 4 ps 1.2 lc 2 title '',\


set out
