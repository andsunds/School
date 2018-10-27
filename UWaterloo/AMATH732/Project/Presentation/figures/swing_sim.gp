
#input header '\small'

#set format xy "$%g$"

##### höjd,bredd
#set size 1,1.4
# Förklaringsrutan, man kan behöva ställa in width för att den ska bli lagom bred
set key samplen 1.5 box top left height .5
#set grid #Sätter på rutnät


#set decimalsign ','#Bestämmer decimalseparator




set terminal epslatex font ',11' color size 16cm,8cm
set key width +.5
set output "swing_sim_energy.tex"
set xlabel '$t$'
set ylabel '$E(t)$'
set logscale y
set format y "$10^{%L}$"

plot [0:600] [1e-3:1.2e-1] "swing_sim_optimal.tsv" using 1:3 with lines lt 1 lw 3 \
lc rgb "black" title '$\phi(t)$',\
"swing_sim_optimal_lin.tsv" using 1:3 with lines lt 1 lw 3 \
lc rgb "blue" title '$\phi_\text{lin}(t)$',\
((0.1**2/8)*exp(3*0.01*x/4)) with lines lt 2 lw 6\
lc rgb "red" title '$E_0\exp(3\epsilon{t}/4)$'



set format y "$%g$"
set ytics -0.8,0.4,0.8
set ytics -1,1,1
set key at 10,0.96 width -1 samplen 1

unset logscale y

#################### plottar ####################

set terminal epslatex font ',12' color size 17cm,8cm
set output "swing_sim.tex"
set xlabel '$t$'
set ylabel 'Original EOM'

set multiplot layout 1,2
#### x-axel     y-axel
plot [0:600] [-1:1] "swing_sim_optimal.tsv" using 1:2 with lines lt 1 lw 2 \
lc rgb "black" title '$\phi(t)$',\
(0.1*exp(3*0.01*x/8)) with lines lt 2 lw 5\
lc rgb "red" title '$\pm a_0\exp(3\epsilon{t}/8)$',\
(-0.1*exp(3*0.01*x/8)) with lines lt 2 lw 5\
lc rgb "red" notitle

set ylabel 'Linearized EOM'
plot [0:600] [-1:1] "swing_sim_optimal_lin.tsv" using 1:2 with lines lt 1 lw 2 \
lc rgb "blue" title '$\phi_\text{lin}(t)$',\
(0.1*exp(3*0.01*x/8)) with lines lt 2 lw 5\
lc rgb "red" title '$\pm a_0\exp(3\epsilon{t}/8)$',\
(-0.1*exp(3*0.01*x/8)) with lines lt 2 lw 5\
lc rgb "red" notitle


