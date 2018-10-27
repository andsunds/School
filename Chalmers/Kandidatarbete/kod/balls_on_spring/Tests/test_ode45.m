%% Test av 2:a ordn ode
clc;clf;clear

x0=[0 10];
k=2;
n=10;
m=1;

odefun=@(x) [x(2); -k*x(1)/m];

[T,X]=ode45(@(t, x) odefun(x), [0,25], x0);

plot(T,X(:,1))

