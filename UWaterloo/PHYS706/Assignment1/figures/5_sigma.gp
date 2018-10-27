set terminal epslatex font ',10' color size 14cm,7cm



#set xzeroaxis
#unset xtics
#set yzeroaxis
#unset ytics

set samples 2000


#set format y ""
#set format x ""


#################### plottar ####################
set output "5_sigma.tex"
set xlabel '$\xi$'
set ylabel '$\Sigma_x(\xi, \theta)$'

#set multiplot layout 2,2

sigma(x,t)=1/((1+x)**2+t**2)-1/((1-x)**2+t**2)


set key samplen 3 box top right height .5
set key width 2
plot [0:4] [-5:5] sigma(x,0.5) lt 1 lw 5 lc 0 title '$\theta=1/2$',\
sigma(x,1) lt 0 lw 5 title '$\theta=1$',\
sigma(x,2) lt 5 lw 5 lc 0 title '$\theta=2$',\
