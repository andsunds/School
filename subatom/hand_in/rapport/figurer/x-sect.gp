set terminal epslatex size 12cm,7cm font ',11' color
#input header '\small'

#set format xy "$%g$"

##### höjd,bredd
#set size 1,1.4
# Förklaringsrutan, man kan behöva ställa in width för att den ska bli lagom bred
set key samplen 2 box b r height .8 horizontal
set grid #Sätter på rutnät
set logscale y #Ställer in rutnät
#set ytics 1e-10,10,1e1
set format y "$10^{%L}$"; set mytics 10#minor ticks
set xlabel '$E$ /[MeV]'
set ylabel '$\sigma$ /[mb]'
set decimalsign ','#Bestämmer decimalseparator

#set key at 9,3.2e-2

#################### plottar ####################
#energydepleted
set key width -3
set output "Se_76_x-sect.tex"

#mapp='../../'
#fil='rp035076.tot'#'../../x-sect/Se_76/rp035076.tot'
s=''#'\footnotesize '

#### x-axel     y-axel
plot [5:10][1e-2:10e2] \
'rp033073.tot' using 1:2 with linespoints title s.'$^{73}_{33}\mathrm{As}$',\
'rp034076.tot' using 1:2 with linespoints pt 6 title s.'$^{76}_{34}\mathrm{Se}$',\
'rp035076.tot' every ::58::99 using 1:2 with linespoints title s.'$^{76}_{35}\mathrm{Br}$',\
'rp035077.tot' using 1:2 with linespoints title s.'$^{77}_{35}\mathrm{Br}$',\
