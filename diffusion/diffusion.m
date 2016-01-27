%% parabolisk diffusion
clc;clear all;clf

N=1000;
L=10;
dx=L/(N+1);
x=linspace(-L/2,L/2, N+1);

n=100;
f0=zeros(length(x), 1);
f0(floor(N/2+1), 1)=n;%'delta-funktion'

%size(f0)

tau=1;
alpha=1;



[T, F] = ode45(@(t,f) par_odefun(f, tau, alpha, dx), [0,10], f0);

figure(1);clf
plot(T, sum(F,2));
figure(2);clf
plot(T, sum(repmat(x.^2,length(T),1).*F(:,1:(N+1)),2)); %variansen

figure(3)
p=plot(x,F(1,:));

playbackspeed=1;
dt=[.2; diff(T)./playbackspeed];

for i=1:length(T)
    pause(dt(i))  
    set(p, 'YData',F(i,:));
end

fprintf('Färdig\n')







%% hyperbolisk diffusion
clc;clear all;clf

N=100;
L=1;
dx=L/(N+1);
x=linspace(-L/2,L/2, N+1);

n=100;
f0=zeros(length(x), 2);
f0(floor(N/2+1), 1)=n;%'delta-funktion'

%size(f0)

tau=10;
alpha=1;



[T, F] = ode45(@(t,f) hyp_odefun(f, tau, alpha, dx), [0,10], reshape(f0,1, []));

figure(1);clf
plot(T, sum(F(:,1:(N+1)),2));%antalet pariklar, denna borde vara konstant
figure(2);clf
plot(T, sum(repmat(x.^2,length(T),1).*F(:,1:(N+1)),2)); %variansen

%
figure(3);clf
p=plot(x,F(1,1:(N+1)));

playbackspeed=1;
dt=[.2; diff(T)./playbackspeed];

for i=1:length(T)
    pause(dt(i))  
    set(p, 'YData',F(i,1:(N+1)));
end

fprintf('Färdig\n\n')