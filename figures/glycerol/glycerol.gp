set terminal epslatex size 8cm,7.5cm font ',10' color
#input header '\small'


set key samplen 2 box height .5 bottom left
set ylabel '$\mathcal{R}$ /[$\unit{s^{-1}}$]'
set xlabel '$\tau^2$ /[$\unit{ms^2}$]'
set grid
set format y '%g'

g_p1=system("head -1 glycerol_fit.tsv | awk '{print $1}'")+0
g_p2=system("head -1 glycerol_fit.tsv | awk '{print $2}'")+0
#t_p1=system("head -1 tap-water_fit.tsv | awk '{print $1}'")+0
#t_p2=system("head -1 tap-water_fit.tsv | awk '{print $2}'")+0

s='\footnotesize '

set grid

set key width -7
set output "glyc.tex"
plot [0:50] [30:60] "glyc-paper.tsv" using 1:(1e3*$2):(1e3*($2-$3)):(1e3*($2+$4)) with yerrorbars lt 1 pt 4 ps 1.2 lc 1 title s.'Gylcerol paper',\
                   "glyc-paper-ex.tsv" using 1:(1e3*$2):(1e3*($2-$3)):(1e3*($2+$4)) with yerrorbars lt 1 pt 8 ps 1.2 lc 2 title s.'Ex. glycerol paper',\
                   "glycerol.tsv" using 1:(1e3*$2):(1e3*($2-$3)):(1e3*($2+$4)) with yerrorbars lt 1 pt 7 ps 1.5 lc 7 title s.'Liquid glycerol',\
                   1e3*(x*g_p1+g_p2) lt 0 lw 10 lc 7 title '',\
                   #1e3*(x*t_p1+t_p2) lt 0 lw 4 lc rgb 'red' title ''




set out
