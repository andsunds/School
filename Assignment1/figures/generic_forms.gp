set terminal epslatex font ',12' color size 16cm,15cm

unset key

set xzeroaxis
unset xtics
set yzeroaxis
unset ytics

set key width 0


set format y ""
set format x ""


#################### plottar ####################
set output "generic_forms.tex"
set xlabel '$M$'
set ylabel '$F(M, t, u)$'

set multiplot layout 2,2

t=1
u=1
set title 'I. $t\ge0$, $u\ge0$'
plot [-1:1][-1:1] (t*x**2+u*x**4+x**6) lt 1 lw 5 lc 0 notitle

t=-2.5
u=2
set title 'II. $t<0$, $u$ any value; or $t=0$, $u<0$'
plot [-1.1:1.1] [-1:1] (t*x**2+u*x**4+x**6) lt 1 lw 5 lc 0 notitle

t=.5
u=-1.6
set title 'III. $0<t<u^2/4$, $u<0$'
plot [-1.3:1.3][-.25:.25](t*x**2+u*x**4+x**6) lt 1 lw 5 lc 0 notitle


t1=.4
u1=-1.2
t2=.6
u2=-1.2
set title 'IV. $t>u^2/4$, $u<0$'

set key samplen 2 box bottom right height .5
set key width -4
plot [-1:1] [-.2:.2] (t1*x**2+u1*x**4+x**6)\
lt 1 lw 5 lc 0 title '$t>u^2/4$, but not by much',\
(t2*x**2+u2*x**4+x**6)\
lt 2 lw 5 lc 0 title '$t\gg u^2/4$',\
