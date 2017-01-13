%% A1, Q2 (b)
clc;clear;clf

x=linspace(-1,1,500);
t1=.1;
t2=.3;
t3=3;
s=1;
a=1;

f=cos(pi*x)+2*cos(3*pi*x);
ft1=cos(pi*(x-t1*a)).*exp(-s*1^2*t1)+2*cos(3*pi*(x-a*t1))*exp(-s*3^2*t1);
ft2=cos(pi*(x-t2*a)).*exp(-s*1^2*t2)+2*cos(3*pi*(x-a*t2))*exp(-s*3^2*t2);
ft3=cos(pi*(x-t3*a)).*exp(-s*1^2*t3)+2*cos(3*pi*(x-a*t3))*exp(-s*3^2*t3);

plot(x,f,'-k'),hold on
plot(x,ft1,'--k')
plot(x,ft2,'-.k')
plot(x,ft3, ':k')

set(gca,'fontsize', 35, 'yTick', [-3,0,3])

ylabel('$u(x, t)$','interpreter','LaTeX','fontsize',35)
xlabel('$x$','interpreter','LaTeX','fontsize',35)

l=legend('$t=0$,', '$t=0.1$,', '$t=0.3$,', '$t=3$');
set(l,'interpreter','latex', 'fontsize',35,...
    'location','south', 'orientation','horizontal', 'linewidth',2)


%% A1, Q2 (a)
clc;clear

A=repmat(-2:2,5,1).^repmat((0:4).',1,5);

b2=[0,0,1,0,0].'*24;%this should be 12...
b3=[0,0,0,1,0].'*12;

A\b2

A\b3

%% A1, Q2 (b)
clc;clear

A=repmat(0:5,6,1).^repmat((0:5).',1,6);
b=[zeros(1,4),1,0].'*24

A\b

%% A1, Q3 (a+b)
clc;clf;clear

f=@(x) sin(x);
Df=@(x) cos(x);
D4f=@(x) sin(x);

FDf=@(dx, x0) (f(x0+dx)-f(x0))./dx;

FD42f=@(dx,x0) (3*f(x0+0*dx)-14*f(x0+1*dx)+26*f(x0+2*dx)-24*f(x0+3*dx)+11*f(x0+4*dx)-2*f(x0+5*dx))./(dx.^4);



dx=logspace(-15,-1,100);
x0=1;

e=abs(Df(x0)-FDf(dx,x0));
e42=abs(D4f(x0)-FD42f(dx,x0));

plot(dx,e,'-'), hold on
plot(dx,e42,'--')
set(gca,'yscale','log', 'xscale','log', 'fontsize',23,...
    'XTick',10.^(-15:-1), 'xLim',[1e-15,1e-1], 'yLim', [1e-10,1e0])

xlabel('$\Delta{x}$','interpreter','LaTeX')
ylabel('Difference between exact and approximated derivatives','interpreter','LaTeX')

l=legend('FD, 1st derivative', 'FD, 4th derivative (from problem 2b)');
set(l, 'location', 'northWest', 'interpreter','LaTeX', 'fontsize', 25)

