%%
clc;clf;clear

eps=0.01;
Lambda=[1,3];
xmin=zeros(size(Lambda));
gmin=zeros(size(Lambda));

linetype={'-b','--r'};

for i=1:length(Lambda)
lambda=Lambda(i);
f=@(x) (x-lambda).^2.*(1-x.^2) +eps*x;
g=@(x) x.^(-2)+eps*(x-lambda).^(-2);

xmin(i)=lambda/(1+eps^(1/3));
gmin(i)=(1+eps^(1/3))^3/lambda^2;

x1=-2;x2=6;
x=linspace(x1,x2,500);
plot(x,g(x),linetype{i})
hold on
end
plot([x1,x2],[1,1],':k',xmin,gmin,'ko','markersize',10)
axis([x1,x2,0,4])

ylabel('$D(x)$','interpreter','latex')
xlabel('$x$','interpreter','latex')

l=legend('$\lambda=1$','$\lambda=3$');
set(l,'interpreter','latex')
set(gca,'fontsize',15)


%%
clc;clear
eps=0.01;
lambda=1-0.00001;

f=@(x) x.^(-2)+eps*(x-lambda).^(-2)-1;

d=eps^(1/2);
x0=lambda+d*lambda/sqrt(lambda^2-1);
abs(f(x0))
d^2

lambda=1;
f=@(x) x.^(-2)+eps*(x-lambda).^(-2)-1;
d=eps^(1/3);
x1=lambda+d*exp(1i*2*pi/3)/2^(1/3);
abs(f(x1))
d^2


%%
n=0:2
abs(f(1+eps^(1/3)*exp(1i*2*pi*n/3)/2^(1/3)))
abs(g(1+eps^(1/3)*exp(1i*2*pi*n/3)/2^(1/3)))
eps^(2/3)

%%
clc;clf

a=cos(2*pi/3)/2^(1/3)
b=sin(2*pi/3)/2^(1/3)
l=linspace(0,0.95);
y=l./sqrt(1-l.^2);
plot(l,y,[0,1],b*[1,1],':k')














