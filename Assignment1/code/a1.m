%% A1, Q2 a)
clc;clear

A=repmat(-2:2,5,1).^repmat((0:4).',1,5);

b2=[0,0,1,0,0].'*24;%this should be 12...
b3=[0,0,0,1,0].'*12;

A\b2

A\b3

%% A1, Q2 a)
clc;clear

A=repmat(0:5,6,1).^repmat((0:5).',1,6);
b=[zeros(1,4),1,0].'*24

A\b

%% A1, Q3 a)
clc;clf;clear

f=@(x) sin(x);
Df=@(x) cos(x);

FDf=@(dx, x0) (f(x0+dx)-f(x0))./dx;

dx=logspace(-15,-1,100);
x0=1;

e=abs(Df(x0)-FDf(dx,x0));

plot(dx,e)
set(gca,'yscale','log', 'xscale','log', 'fontsize',20,...
    'XTick',10.^(-15:-1), 'xLim',[1e-15,1e-1])

xlabel('$\Delta{x}$','interpreter','LaTeX')
ylabel('$|FDf(\Delta{x}, x_0)-Df(x)|$','interpreter','LaTeX')

