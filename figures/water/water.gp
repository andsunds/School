set terminal epslatex size 8cm,7.5cm font ',10' color
#input header '\small'


set key samplen 2 box height .5 bottom right
set ylabel '$\mathcal{R}$ /[$\unit{s^{-1}}$]'
set xlabel '$\tau^2$ /[$\unit{ms^2}$]'
set grid
set format y '%g'

d_p1=system("head -1 d-water_fit.tsv | awk '{print $1}'")+0
d_p2=system("head -1 d-water_fit.tsv | awk '{print $2}'")+0
t_p1=system("head -1 tap-water_fit.tsv | awk '{print $1}'")+0
t_p2=system("head -1 tap-water_fit.tsv | awk '{print $2}'")+0

s='\footnotesize '

set grid

set key width -3.5
set output "water.tex"
plot [0:2500] [0:10] "d-water.tsv" using 1:(1e3*$2):(1e3*($2-$3)):(1e3*($2+$4)) with yerrorbars lt 1 pt 3 ps 1.2 lc 7 title s.'On mark',\
                       "tap-water.tsv" using 1:(1e3*$2):(1e3*($2-$3)):(1e3*($2+$4)) with yerrorbars lt 1 pt 9 ps 1.2 lc 6 title s.'Off mark',\
                       1e3*(x*d_p1+d_p2) lt 0 lw 4 lc 7 title '',\
                       1e3*(x*t_p1+t_p2) lt 0 lw 4 lc 6 title ''




set out
