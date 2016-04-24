%% del 1, mäta intensitet vid de olika positinoerna
clc;clearvars


W1=1;
W2=(2.5/5.5)^2;
W3=(2.5/8.5)^2;


b=14.1/60; %Bq

W=[W3; W2; W3; W2; W1];

N=[899; 429; 753; 1984; 2970];

Int=(N*2-b)./W;

disp(Int.'/1e3)

%% del 2, Halveringstid
clc, clf;clearvars

N=[1120,714, 531,401,320,265,237,210,176,140,146,278,207,209,134,143,147,108,105,93,80,73,60,109,94,64];
dt=[10*ones(1,11), 20*ones(1,12), 50*ones(1,3)]/60;
t=(cumsum([0,10+dt(1:end-1)*60])+30)/60+.5*dt;
b=14.1;
A=N./dt-b;

plot(t,A, '-o'), hold on


I2=(15:26);
C2=[ones(size(I2.')), -t(I2).']\log(A(I2).');

T_halv_110= log(2)/C2(2)*60





I0=1:8;
plot(t(I0), A(I0)-exp(C2(1)-C2(2)*t(I0)), 'o-')

I1=(1:7);
C1=[ones(size(I1.')), -t(I1).']\log((A(I1)-exp(C2(1)-C2(2)*t(I1))).');

T_halv_108= log(2)/C1(2)*60

%plottar anpassnigar
x=linspace(0,14);
plot(x, exp(C2(1)-C2(2)*x))
plot(x, exp(C1(1)-C1(2)*x))

axis([0,14, 1e1, 1e4])

set(gca,'fontsize', 20, 'yscale', 'log', 'xscale', 'lin')

%% överlagrade exponentialer
clc;clf;clearvars

x=linspace(-10,10);
k1=1;%1.1540;
l1=.5;%;1.7414;
k2=1e-0;%0.2624;
l2=.1;%0.2856;


plot(x, k1*exp(-l1*x), x, k2*exp(-l2*x)), hold on
plot(x, k1*exp(-l1*x)+k2*exp(-l2*x))
plot(x, (k1*exp(l1*x)+k2*exp(l2*x)).^-1)

axis([-8,8, 1e-2, 1e2])

set(gca,'fontsize', 20, 'yscale', 'log', 'xscale', 'lin')










