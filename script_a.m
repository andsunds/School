%%
clc;clear all
p='Data/';addpath(p);

files_t=dir([p,'Pb','*_t.lvm']);
files_v=dir([p,'Pb','*_v.lvm']);

l=length(files_t);

T=zeros(1,l);
Pb_order = [2 3 4 6 5]-1;
Pb_x = [0,1,2,3,4].';
beta1 = zeros(l,1);
gamma1 = zeros(l,1);

for i=1:l
    time = load(files_t(i).name);
    time = (time-min(min(time)))/1000;
    volt = load(files_v(i).name);
    cels = V2C(volt,Pb_order);
    
    
    T(i)=str2double(files_t(i).name(regexp(files_t(i).name, '\d')));
    
    
    
    
    
    [log_r1, phi1] = Heterodyn(time,cels,T(i)*60, 4, Pb_order);
    [log_r2, phi2] = Heterodyn(time,volt,T(i)*60, 4, Pb_order);

    
    X1 = [ones(3,1),Pb_x(1:3)]\[log_r1(1:3),phi1(1:3)];
    log_V0 = X1(1,1);
    phi0 = X1(1,2);
    gamma1(i) = -X1(2,1);
    beta1(i) = -X1(2,2);
    
    N=5;
    X2 = [ones(N,1),Pb_x(1:N)]\[log_r2(1:N),phi2(1:N)];
    gamma2(i) = -X2(2,1);
    beta2(i) = -X2(2,2);
    
    
%     subplot(2,2,[1,3])
%     plot(repmat(time(:,1)/60,1,5),cels)
%     title(sprintf('Periodtid: %dmin',T(i)))
%     xlabel('Tid [min]')
%     ylabel('Temperatur [C]')
%     
%     subplot(2,2,2)
%     plot(Pb_x,log_r1,'o')
%     hold on
%     plot(Pb_x, log_V0 - Pb_x*gamma1(i),'-')
%     hold off
%     subplot(2,2,4)
%     plot(Pb_x,phi1,'v')
%     hold on
%     plot(Pb_x, phi0 - Pb_x*beta1(i),'-')
%     hold off
%     pause(.1)
end


clf
subplot(1,2,1)
tau = T.^(-1/2);
plot(tau,beta1,'bo'),hold on
plot(tau,beta2,'rx')
grid on
%plot(tau,gamma,'rv')
hold off
%axis([0,max(tau)*1.1,0,max(beta1)*1.1])
title('beta')
legend('temperatur', 'spänning', 'location', 'NorthWest')


subplot(1,2,2)
plot(tau,gamma1,'bv'), hold on
plot(tau,gamma2,'r^')
grid on
%axis([0,max(tau)*1.1,0,max(gamma1)*1.1])
title('gamma')
legend('temperatur', 'spänning', 'location', 'NorthWest')





























