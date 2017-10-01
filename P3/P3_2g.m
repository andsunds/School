%%
clc;clf;clear
eps=0.001;
lambda=1;


f=@(x) (x-lambda).^2.*(1-x.^2) +eps*x;
g=@(x) f(x)./(x-lambda).^2.*(1-x.^2);
n=0:2;

abs(f(1+eps^(1/3)*exp(1i*2*pi*n/3)/2^(1/3)))
abs(g(1+eps^(1/3)*exp(1i*2*pi*n/3)/2^(1/3)))
eps^(2/3)

%%
clc

cos(2*pi/3)/2^(1/3)
sin(2*pi/3)/2^(1/3)