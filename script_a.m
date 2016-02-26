%%
clc; clf;clear all
p='Data/';addpath(p);

files_t=dir([p,'Pb','*_t.lvm']);
files_v=dir([p,'Pb','*_v.lvm']);

l=length(files_t);

T=zeros(l,1);
Pb_order = [2 3 4 6 5]-1;
Pb_x = [0,1,2,3,4].';
beta1  = zeros(l,1); beta2  = zeros(l,1);
gamma1 = zeros(l,1); gamma2 = zeros(l,1);

for i=1:l
    %j=floor((1+i)/2);%övertoner
    j=i;
    name_t=files_t(j).name;
    name_v=files_v(j).name;
    time = load(name_t);
    time = (time-min(min(time)))/1000;
    volt = load(name_v);
    cels = V2C(volt,Pb_order);
    T(i)=str2double(name_t(regexp(name_t, '\d')));
    
    %cels = highpass(time,cels,0.1/T(j));
    
    Nper=5;
    
    if i==j*2
        T(i)=T(i)/3;%tredje övertonen
        Nper=Nper*3;
        
    end
    
    
    
    [log_r1, phi1] = Heterodyn(time,cels,T(i)*60, Nper, Pb_order);
    [log_r2, phi2] = Heterodyn(time,volt,T(i)*60, Nper, Pb_order);
    
    N=3;
    
    %%%%%%Roetzels metod%%%%%%
    %(phi1(:)-phi1(1))
    %k=-(phi1(1:N)-phi1(1))\(log_r1(1:N)-log_r1(1))
    %k=mean((log_r1(2:4)-log_r1(1))./(phi1(2:4)-phi1(1)))
    %tau=((1-k.^2)./(2*k))*T(i)/2/pi
    
    
    X1 = [ones(N,1),Pb_x(1:N)]\[log_r1(1:N),phi1(1:N)];
    log_V0 = X1(1,1);
    phi0 = X1(1,2);
    gamma1(i) = -X1(2,1);
    beta1(i) = -X1(2,2);
    
    
    X2 = [ones(N,1),Pb_x(1:N)]\[log_r2(1:N),phi2(1:N)];
    gamma2(i) = -X2(2,1);
    beta2(i) = -X2(2,2);
    
    
    subplot(2,2,[1,3])
    plot(repmat(time(:,1)/60,1,5),cels)
    title(sprintf('Periodtid: %dmin',T(i)))
    xlabel('Tid [min]')
    ylabel('Temperatur [C]')
    
    subplot(2,2,2)
    plot(Pb_x,log_r1,'o')
    hold on
    plot(Pb_x, log_V0 - Pb_x*gamma1(i),'-')
    hold off
    subplot(2,2,4)
    plot(Pb_x,phi1,'v')
    hold on
    plot(Pb_x, phi0 - Pb_x*beta1(i),'-')
    hold off
    pause(.2)
end
%%
clf
w = 2*pi*T.^(-1);

a=3.4e-1;
tau=.01;
x=linspace(0,4);

BETA=@(tau, a) sqrt(w/(2*a)).*sqrt(sqrt((tau*w).^2+1)+tau*w);

GAMMA=@(tau,a) sqrt(w/(2*a)).*sqrt(sqrt((tau*w).^2+1)-tau*w);

X=fminsearch(@(X) sum( (BETA(X(1), X(2)) - beta1).^2 ), [1, 1])
Y=fminsearch(@(X) sum( (GAMMA(X(1), X(2)) - gamma1).^2 ), [1, 1])

tau=X(1);a=X(2);
y1=sqrt(x/(2*a)).*sqrt(sqrt((tau*x).^2+1)+tau*x);
tau=Y(1);a=Y(2);
y2=sqrt(x/(2*a)).*sqrt(sqrt((tau*x).^2+1)-tau*x);


subplot(1,2,1)
plot(w,beta1,'bo'),hold on
plot(w,beta2,'rx')
plot(x,y1)
grid on
hold off

title('beta')
legend('temperatur', 'spänning', 'location', 'Best')





subplot(1,2,2)
plot(w,gamma1,'bv'), hold on
plot(w,gamma2,'r^')
plot(x,y2)
grid on
hold off
%axis([0,max(tau)*1.1,0,max(gamma1)*1.1])
title('gamma')
legend('temperatur', 'spänning', 'location', 'best')

%%
pause()

clc;clf

plot(w,beta1./gamma1,'bo'),hold on
xlabel('omega [s^{-1}]')
ylabel('beta./gamma')
grid on


%w=tau.^2;
%y= @(t, w)(sqrt(sqrt((w*t).^2+1)-t*w)./sqrt(sqrt((w*t).^2+1)+t*w));

%t0=fminbnd(@(t)sum((y(t, w)-(beta1./gamma1).').^2), 0,1000)
%W=linspace(0.02, 0.2);
%plot(W, y(t0, W))


%% Test av h�gpassfilter
clc; clf;clear all
p='Data/';addpath(p);

files_t=dir([p,'Pb','*_t.lvm']);
files_v=dir([p,'Pb','*_v.lvm']);

l=length(files_t);

T=zeros(1,l);
Pb_order = [2 3 4 6 5]-1;
Pb_x = [0,1,2,3,4].';

for i=1:l
    %j=floor((1+i)/2);%övertoner
    j=i;
    name_t=files_t(j).name;
    name_v=files_v(j).name;
    time = load(name_t);
    time = (time-min(min(time)))/1000; %Konvertera till sekunder
    volt = load(name_v);
    cels = V2C(volt,Pb_order);
    T(j) = 60*sscanf(name_t,'Pb_%dmin'); % Periodtid i minuter
    
    cels_highpass = highpass(time,cels,0.1/T(j));
    
    subplot(2,1,1)
    plot(repmat(time(:,1)/60,1,5),cels)
    title(sprintf('Periodtid: %dmin',T(i)/60))
    xlabel('Tid [min]')
    ylabel('Temperatur [C]')
    
    subplot(2,1,2)
    plot(repmat(time(:,1)/60,1,5),cels_highpass)
    title('H�gpassfiltrerad signal')
    xlabel('Tid [min]')
    ylabel('Temperatur [C]')
    
    pause()
end

