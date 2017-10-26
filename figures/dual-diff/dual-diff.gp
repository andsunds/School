set terminal epslatex size 15cm,7.5cm font ',10' color
#input header '\small'


set key samplen 2 box height .5 top right
set ylabel '$\mathcal{R}$ /[$\unit{s^{-1}}$]'
set xlabel '$\tau^2$ /[$\unit{ms^2}$]'
set grid
set format y '%g'

#d_p1=system("head -1 d-water_fit.tsv | awk '{print $1}'")+0
#d_p2=system("head -1 d-water_fit.tsv | awk '{print $2}'")+0


s='\footnotesize '

set grid

set key width -3.6
set output "diff.tex"
plot [0:50] [0:60] "wet-pp-ex.tsv" using 1:(1e3*$2):(1e3*($2-$4)):(1e3*($2+$5)) with yerrorbars lt 1 pt 4 ps 1.2 lc 7 title s.'Ex. wet printer paper',\
                   "wet-pp-ex.tsv" using 1:(1e3*$3):(1e3*($3-$6)):(1e3*($3+$7)) with yerrorbars lt 1 pt 5 ps 1.2 lc 7 title '',\
                   "wet-tp.tsv" using 1:(1e3*$2) lt 1 pt 6 ps 1.2 lc 6 title s.'Wet tissue paper',\
                   "wet-tp.tsv" using 1:(1e3*$3) lt 1 pt 7 ps 1.2 lc 6 title '',\
                   #"wet-pp-ex.tsv" using 1:(1e3*$4) lt 1 pt 12 ps 1.2 lc 7 title '',\
                   #"wet-tp.tsv" using 1:(1e3*$4) lt 1 pt 21 ps 1.2 lc 6 title '',\
                       #1e3*(x*d_p1+d_p2) lt 0 lw 4 lc rgb 'blue' title '',\
                       #1e3*(x*t_p1+t_p2) lt 0 lw 4 lc rgb 'red' title ''
#Need to add fast and slow mode



set out
