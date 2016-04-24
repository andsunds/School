set terminal epslatex font ',11' size 12cm,8cm

#set format xy "$%g$"

##### höjd,bredd
set size 1,1
# Förklaringsrutan, man kan behöva ställa in width för att den ska bli lagom bred
set key samplen 2 box outside#height .4

set grid #Sätter på rutnät

set xlabel '$\varLambda$'
set ylabel 'MSD $/\left[\text{godt. längdenhet}^2\right]$'
set decimalsign ','#Bestämmer decimalseparator


#################### plottar ####################
#kvot
set key width -4
set output "isotropi_kvot.tex"
set logscale xy #skala på axlarna
set format y '$10^{%L}$'

fil="egenvardeskvot.tsv"

#### x-axel     y-axel
plot [1:100] [1e-2:1] fil using 1:2 with lines \
lt -1 lw 2 lc rgb "gray" title 'kvot',\
fil using 3:4 with lines \
lt 3 lw 2 lc rgb "blue" title 'brus 14',\
fil using 5:6 with lines \
lt 7 lw 2 lc rgb "blue" title 'brus 9',\
fil using 7:8 with linespoints \
lt 1 lw 2 lc rgb "black" title 'dvala',\
fil using 9:10 with linespoints \
lt 10 lw 2 lc rgb "black" title 'aktiv',\
fil using 11:12 with lines \
lt 5 lw 2 lc rgb "red" title 'spring 2,6\,s',\
fil using 13:14 with lines \
lt 8 lw 2 lc rgb "red" title 'spring 6,0\,s',\


#asymmetri
set output "isotropi_asymmetri.tex"
set key top right
unset logscale xy
set format y "$%g$"
set ytics 0,0.2,1
set xtics 0,0.2,1


fil="asymmetri.tsv"
#### x-axel     y-axel
plot [0:1] [0:1] fil using 1:2 with lines \
lt -1 lw 2 lc rgb "gray" title 'kvot',\
fil using 3:4 with lines \
lt 3 lw 2 lc rgb "blue" title 'brus 14',\
fil using 5:6 with lines \
lt 7 lw 2 lc rgb "blue" title 'brus 9',\
fil using 7:8 with linespoints \
lt 1 lw 2 lc rgb "black" title 'dvala',\
fil using 9:10 with linespoints \
lt 10 lw 2 lc rgb "black" title 'aktiv',\
fil using 11:12 with lines \
lt 5 lw 2 lc rgb "red" title 'spring 2,6\,s',\
fil using 13:14 with lines \
lt 8 lw 2 lc rgb "red" title 'spring 6,0\,s',\
