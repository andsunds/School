set terminal epslatex
#energydepleted
set output "MSD_ed.tex"
#set format xy "$%g$"
#Comment
set xlabel '$\Delta{t}$ /[s]'
set ylabel 'MSD /[godt. längdenhet$^2$]'
set key samplen 2 box width -3 top left height .4
set decimalsign ','
set grid
set logscale xy
set ytics 1e-9,10,1e-7;set format y "$10^{%L}$"


plot [1e-2:1e1] [1e-9:2e-7] "MSD.tsv" using 1:2 with lines lt 1 lc rgb "blue" title '$S(\Delta{t})$ energydepleted',\
"MSD.tsv" using 1:3 with lines lt 2 lc rgb "red" title '$s(\Delta{t})$ energydepleted'

#logphase
set output "MSD_lp.tex"
set key width -1.5

plot [1e-2:1e1] [1e-9:2e-7] "MSD.tsv" using 1:4 with lines lt 1 lc rgb "blue" title '$S(\Delta{t})$ logphase',\
"MSD.tsv" using 1:5 with lines lt 2 lc rgb "red" title '$s(\Delta{t})$ logphase'
