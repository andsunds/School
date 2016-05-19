set terminal epslatex size 12cm,7cm font ',10' color
#input header '\small'

#set format xy "$%g$"

##### höjd,bredd
#set size 1,1.4
# Förklaringsrutan, man kan behöva ställa in width för att den ska bli lagom bred
set key samplen 2 box t l height .8
set grid #Sätter på rutnät
#set logscale y #Ställer in rutnät
#set ytics 1e-10,10,1e1
#set format y "$10^{%L}$"; set mytics 10#minor ticks
set xlabel '$E$ /[MeV]'
set ylabel '$\sigma$ /[b]'
set decimalsign ','#Bestämmer decimalseparator

#set key at 9,3.2e-2

#################### plottar ####################
#energydepleted
set key width -3
set output "Se_76_x-sect.tex"

fil='rp035076.tot'#'../../x-sect/Se_76/rp035076.tot'
s=''#'\footnotesize '

#### x-axel     y-axel
plot [5:10][0:700] \
fil using 1:2 with linespoints lt 1 pt 5 \
lw 3 lc rgb "blue" title s.'$^{76}_{34}\mathrm{Se}$'



