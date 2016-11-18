
#input header '\small'

#set format xy "$%g$"

##### höjd,bredd
#set size 1,1.4
# Förklaringsrutan, man kan behöva ställa in width för att den ska bli lagom bred
set key samplen 1.5 box top left height .5
#set grid #Sätter på rutnät
#set logscale y

set format y "$%g$"
set ytics -0.8,0.4,0.8
set ytics -1,1,1

#set decimalsign ','#Bestämmer decimalseparator




#################### plottar ####################
#energydepleted std
set terminal epslatex font ',12' color size 15cm,6cm
set key width -2
set output "swing_sim.tex"
set xlabel '$t$'
set ylabel 'Original EOM'
s=''
set multiplot layout 1,2
#### x-axel     y-axel
plot [0:600] [-1:1] "swing_sim_optimal.tsv" using 1:2 with lines lt 1 \
lc rgb "black" title s.'$\phi(t)$',\
(0.1*exp(3*0.01*x/8)) with lines lt 2 lw 5\
lc rgb "red" title s.'$\pm a_0\exp(3\epsilon{t}/8)$',\
(-0.1*exp(3*0.01*x/8)) with lines lt 2 lw 5\
lc rgb "red" notitle,\

set ylabel 'Linearized EOM'
plot [0:600] [-1:1] "swing_sim_optimal_lin.tsv" using 1:2 with lines lt 1 \
lc rgb "blue" title s.'$\phi_\text{lin}(t)$',\
(0.1*exp(3*0.01*x/8)) with lines lt 2 lw 5\
lc rgb "red" title s.'$\pm a_0\exp(3\epsilon{t}/8)$',\
(-0.1*exp(3*0.01*x/8)) with lines lt 2 lw 5\
lc rgb "red" notitle,\



