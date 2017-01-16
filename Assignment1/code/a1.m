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
%N.B! 
%The row numbers, from the declaration of f to e42, should be kept at 55-69!
f=@(x) sin(x); %the original function
Df=@(x) cos(x); %the first derivative (exact)
D4f=@(x) sin(x); %the fourth derivative (exact)

%FD first derivative
FDf=@(dx, x0) (f(x0+dx)-f(x0))./dx; 
%FD fourth derivative
FD42f=@(dx,x0) ( 3*f(x0+0*dx)-14*f(x0+1*dx)+26*f(x0+2*dx)...
                 -24*f(x0+3*dx)+11*f(x0+4*dx)-2*f(x0+5*dx) )./(dx.^4);
%initializing
dx=logspace(-15,-1,100);
x0=1;
%Calculating the absolute errors
e=abs(Df(x0)-FDf(dx,x0));
e42=abs(D4f(x0)-FD42f(dx,x0));

y11=.4*dx.^1;
y12=3e-17*dx.^(-1);
y21=2.5*dx.^2;
y22=3e-15*dx.^(-4);


plot(dx,e,'+-','linewidth',3,'markersize',8), hold on
plot(dx(50:end),e42(50:end),'*-', 'linewidth',3,'markersize',10)

plot(dx,y11,'-k', dx,y12,'--k', 'linewidth',2)
plot(dx,y21,'-.k', dx,y22,':k', 'linewidth',2)

set(gca,'yscale','log', 'xscale','log', 'fontsize',23,...
    'XTick',10.^(-15:-1), 'xLim',[1e-15,1e-1], 'yLim', [1e-10,1e1])

xlabel('$\Delta{x}$','interpreter','LaTeX')
ylabel('Difference between exact and approximated derivatives','interpreter','LaTeX')

l=legend('FD 1st derivative', 'FD 4th derivative (from problem 2b)',...
         '$0.4\,\Delta{x}$','$3\times10^{-17}\,\Delta{x}^{-1}$','$2.5\,\Delta{x}^2$','$3\times10^{-15}\,\Delta{x}$');
set(l, 'location', 'north', 'interpreter','LaTeX', 'fontsize', 25)



%% A1, Q4: num sol of heat equation
clc;clf;clear

sigma=1;
dX=0.01; %mesh size
r=.4;
dt=(r*dX^2/sigma)

X=-1:dX:1; %the mesh
N=length(X);

T=.2;
n=floor(T/dt);

u0=sin(pi*X).'; %IC
%u_end=0;      %BC


%The transfer matrix
A=(1-2*r)*eye(N)...%the diagonal
   +[zeros(N-1,2),[zeros(1,N-2); r*eye(N-2)];zeros(1,N)]...%the upper off-diagonal
   +[zeros(1,N);[r*eye(N-2);zeros(1,N-2)],zeros(N-1,2)];%the lower off-diagonal


U=[u0,zeros(N,n)];%init

for i=2:n
    U(:,i)=A*U(:,i-1);
end

%plot(X,U(:,[1,floor(n/3),floor(2*n/3),n]))
plot(X,U(:,1),':'), hold on
plot(X,U(:,floor(n/3)),'-')
plot(X,U(:,floor(2*n/3)),'--')
plot(X,U(:,floor(n)),'-.')


set(gca,'fontsize',30)

xlabel('$x$','interpreter','LaTeX')
ylabel('$U_j^n$','interpreter','LaTeX')

l=legend(sprintf('$t=%0.2f$',0), sprintf('$t=%0.2f$',n/3*dt),...
         sprintf('$t=%0.2f$',2*n/3*dt),sprintf('$t=%0.2f$',n*dt));

set(l, 'location', 'northwest', 'interpreter','LaTeX', 'fontsize', 25)




































