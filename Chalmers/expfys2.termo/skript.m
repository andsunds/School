%% Ladda in data, Aluminium
clear all; clf;%clc;
Al_order = [4,3,5,2,6]-1;
Al_x = [0,5,10,15,20];
addpath('Data')
volt = load('Al_30min_v.lvm', '-ascii'); % volt
time = load('Al_30min_t.lvm', '-ascii'); % millisekunder
% M�tningen skedde vid centrumtiden (size(time)=[:,2])
time = mean(time,2);
% Konvertera till minuter och starta fr�n 0
time = (time-min(min(time)))/60000;
% Konvertera sp�nning till temperatur
cels = V2C(volt,Al_order);
% Anv�nd r�tt ordning p� kanalerna
cels = cels(:,Al_order);

% Ber�kna frekvens och DFT
%dw*T = 2*pi
%time = linspace(0,100,1024);
%cels = sin(time'*2*pi/10);
N = length(time);
f = (0:N-1)/(time(end)-time(1));
y = abs(fft(cels(:,1),size(cels,1)))/37720;
plot(f, y,'-')
%plot(time, cels,'-')
%grid on
axis([0,5/5,0,1])

%legend('0 cm', '5 cm', '10 cm', '15 cm', '20 cm')

tmp = [time,cels,f',y];
%save('Plots/timeseries.tsv', 'tmp', '-ascii','-tabs')

%%

[log_r1, phi1] = Heterodyn(time1,volt1,5*60,Al_order);
subplot(2,1,1)
plot(Al_x,log_r1,'o')
subplot(2,1,2)
plot(Al_x,phi1,'v')

%%
[log_r2, phi2] = Heterodyn(time2,volt2,10*60,Al_order);
subplot(2,1,1)
plot(Al_x,log_r2,'o')
subplot(2,1,2)
plot(Al_x,phi2,'v')



%% Ladda in data, blykulor
clc; clear all;
volt1 = load('Pb_33min_v.lvm', '-ascii');
time1 = load('Pb_33min_t.lvm', '-ascii');
% Tid i sekunder
time1 = (time1-min(min(time1)))/1000;

subplot(2,1,1)
plot(time1(:,1:5),volt1,'-')

volt2 = load('Pb_42min_v.lvm', '-ascii');
time2 = load('Pb_42min_t.lvm', '-ascii');
time2 = (time2-min(min(time2)))/1000;

subplot(2,1,2)
plot(time2(:,1:5),volt2,'-')

Pb_order = [2 3 4 6 5]-1;
Pb_x = [0,1,2,3,4];
%%

[log_r1, phi1] = Heterodyn(time1,volt1,12*60,5,Pb_order);
subplot(2,1,1)
plot(Pb_x,log_r1,'o')
subplot(2,1,2)
plot(Pb_x,phi1,'v')

%%
[log_r2, phi2] = Heterodyn(time2,volt2,21*60,7,Pb_order);
subplot(2,1,1)
plot(Pb_x,log_r2,'o')
subplot(2,1,2)
plot(Pb_x,phi2,'v')

%%
clc; clear all;
%Al_order = [4,3,5,2,6]-1;
%Al_x = [0,5,10,15,20].';
Pb_order = [2 3 4 6 5]-1;
Pb_x = [0,1,2,3,4].';
%T = [2,3,5,7,10,12,17,22,30]; %min
T = [5,7,12,15,21,33,42];
beta = zeros(length(T),1);
gamma = zeros(length(T),1);
for i=1:length(T)
    %filename_t = sprintf('Al_%dmin_t.dat',T(i));
    %filename_v = sprintf('Al_%dmin_v.dat',T(i));
    filename_t = sprintf('Pb_%dmin_t.dat',T(i));
    filename_v = sprintf('Pb_%dmin_v.dat',T(i));
    
    time = load(filename_t,'-ascii');
    time = (time-min(min(time)))/1000;
    volt = load(filename_v,'-ascii');
    
    subplot(2,2,[1,3])
    plot(repmat(time(:,1)/60,1,5),volt)
    title(sprintf('Periodtid: %dmin',T(i)))
    xlabel('Tid [min]')
    ylabel('Spänning [V]')
    
    [log_r, phi] = Heterodyn(time,volt,T(i)*60, 4, Pb_order);
    
    % A*X = B
    % B = [log_r,phi]
    % A = [ones(5,1),Al_x]
    % X = [log_V0,phi0; gamma, beta]
    X = [ones(5,1),Pb_x]\[log_r,phi];
    X = [ones(3,1),Pb_x(1:3)]\[log_r(1:3),phi(1:3)];
    log_V0 = X(1,1);
    phi0 = X(1,2);
    gamma(i) = -X(2,1);
    beta(i) = -X(2,2);
    
    subplot(2,2,2)
    plot(Pb_x,log_r,'o')
    hold on
    plot(Pb_x, log_V0 - Pb_x*gamma(i),'-')
    hold off
    subplot(2,2,4)
    plot(Pb_x,phi,'v')
    hold on
    plot(Pb_x, phi0 - Pb_x*beta(i),'-')
    hold off
    
    pause()
end

%%
subplot(1,2,1)
tau = T.^(-1/2);
plot(tau,beta,'o')
hold on
plot(tau,gamma,'rv')
hold off
axis([0,max(tau)*1.1,0,max(beta)*1.1])
title('beta')

subplot(1,2,2)
plot(tau,gamma,'rv')
axis([0,max(tau)*1.1,0,max(gamma)*1.1])
title('gamma')