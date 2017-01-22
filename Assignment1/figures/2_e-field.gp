set terminal epslatex font ',10' color size 14cm,7cm



#set xzeroaxis
#unset xtics
#set yzeroaxis
#unset ytics

set samples 2000


#set format y ""
#set format x ""


#################### plottar ####################
set output "2_e-field.tex"
set xlabel '$\xi=(r/a)$'
set ylabel '$e(\xi)=E_r(r)\frac{4\pi\enaught a^2}{Q}$'

#set multiplot layout 2,2

cond(x)=x<=1 ? 0 : x**(-2)
n0(x)=x<1 ? x : x**(-2)
nn2(x)=x<1 ? x**(-1) : x**(-2)
np2(x)=x<1 ? x**(3) : x**(-2)


set key samplen 3 box top right height .5
set key width 0
plot [0:2] [0:2] cond(x) lt 1 lw 5 lc 0 title 'Conducting ball ($n\to\infty$)',\
n0(x) lt 0 lw 5 title '$n=0$',\
nn2(x) lt 2 lw 5 lc 0 title '$n=-2$',\
np2(x) lt 5 lw 5 lc 0 title '$n=2$',\
