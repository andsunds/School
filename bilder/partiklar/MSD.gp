set terminal epslatex size 7.5cm,8cm font ',10'\
#input header '\small'

#set format xy "$%g$"

##### höjd,bredd
#set size 1,1.4
# Förklaringsrutan, man kan behöva ställa in width för att den ska bli lagom bred
set key samplen 2 box top left height .5
set grid #Sätter på rutnät
set logscale xy #Ställer in rutnät
set ytics 1e-9,10,1e-7; set format y "$10^{%L}$"
set xlabel '$\Delta{t}$ /[s]'
set ylabel '$1-\text{CDF}$'
set decimalsign ','#Bestämmer decimalseparator


#################### plottar ####################
#energydepleted
set key width -6
set output "MSD_ed.tex"

s='\footnotesize '

#### x-axel     y-axel
plot [1e-2:1e1] [1e-9:2e-7] "MSD.tsv" using 1:2 with lines lt 1 \
lc rgb "blue" title s.'$S(\Delta{t})$ energydepleted',\
"MSD.tsv" using 1:3 with lines lt 5 \
lc rgb "red" title s.'$s(\Delta{t})$ energydepleted',\
1.8e-8*x**0.63 lt 0 title s.'Anpassning $\propto(\Delta{t})^{0,65}$'




#logphase
set output "MSD_lp.tex"
#set key width -3
#### y-axel     x-axel
plot [1e-2:1e1] [1e-9:2e-7] "MSD.tsv" using 1:4 with lines lt 1 \
lc rgb "blue" title s.'$S(\Delta{t})$ logphase',\
"MSD.tsv" using 1:5 with lines lt 5 \
lc rgb "red" title s.'$s(\Delta{t})$ logphase',\
1.8e-8*x**0.8 lt 0 title s.'Anpassning $\propto(\Delta{t})^{0,8}$'
