%PSD

clc;
%clf;
clear all
load('filnamn.mat')
% data_energydepleted =load(filnamn{1});
% data_logphase=load(filnamn{2});

fil=2; %1=ed, 2=log

data=load(filnamn{fil});
C = separera(data);
n=length(C);%antal partiklar

Fs = 100;
n_max=914;
t = C{1}(1:n_max,1);
PSDx=zeros(458,n);
PSDy=zeros(458,n);

load('filnamn.mat')
load(kompl,...
      'intensitet', 'medelsteg', 'sigma_brus', '-mat')

I=intensitet{fil};

%Detta är just nu samma normering som för rörligheten...
lambda=medelsteg{fil}-0e-9; %Räknar ej med bakgrundsbrus
%beräkna koefficienter för lutning i loglog
koef=[ones(size(I)), log(I)]\log(lambda);
koef(1)=exp(koef(1));%konverterar till potenssamband

%figure(3)
%plot(I,medelsteg{fil}-1e-9,'o')
%x=logspace(-1,1.3);
%y=koef(1)*x.^(koef(2));
%hold on 
%plot(x,y)
%set(gca,'xscale','log','yscale','log')


for i=1:n

dx=diff(C{i}(1:n_max+1,2));
dy=diff(C{i}(1:n_max+1,3));
%r=sqrt(dx.^2+dy.^2);
%r=r/(koef(1)*C{i}(1,4)^(koef(2)));
dx=dx/(koef(1)*C{i}(1,4)^(koef(2)));
dy=dy/(koef(1)*C{i}(1,4)^(koef(2)));

N = length(dx);
xdft = fft(dx); %Fouriertransformera
xdft = xdft(1:N/2+1); %Ta bara halva spektrumet
psdx = (1/(Fs*N)) * abs(xdft).^2; %PSD, frekv
%psdx = (1/(2*pi*N)) * abs(xdft).^2; %PSD, vinkelfrekv
psdx(2:end-1) = 2*psdx(2:end-1); %För att bevara totala effekten, vi tog ju bort halva ovan
%0 och Nyquist frek förekommer bara 1 gång därav ej *2

ydft = fft(dy); %Fouriertransformera
ydft = ydft(1:N/2+1); %Ta bara halva spektrumet
psdy = (1/(Fs*N)) * abs(ydft).^2; %PSD
%psdy = (1/(2*pi*N)) * abs(ydft).^2; %PSD, vinkelfrekv
psdy(2:end-1) = 2*psdy(2:end-1);

PSDx(:,i)=psdx;
PSDy(:,i)=psdy;

%freq = 0:Fs/length(x):Fs/2;
%plot(freq,(psdx))
%axis([0 50 0 1e-9])
%pause(1)

%axis([0 50 -210 -150])
end

PSD_sum=(sum(PSDx,2)+sum(PSDy,2))/(2*n_max); %Totala PSD:n


freq = (0:Fs/length(dx):Fs/2)*(2*pi); %Vinkelfrekvens läng x-axeln
%freq = 0:(2*pi)/N:pi; %Vinkelfrekv, 0-pi

%figure(5)
hold on
%plot(freq,10*log10(PSD_sum))
plot(freq,PSD_sum)
grid on
title('PSD för stegen')
%set(gca,'xscale','log','yscale','log')
%xlabel('Frequency (rad/s)')
%ylabel('Power/Frequency (dB/(rad/s))')

mean_PSD=mean(PSDx,2);
std_PSDx=std(PSDx,0,2);
std_PSDy=std(PSDy,0,2);
std_PSD=sqrt(std_PSDx.^2+std_PSDy.^2)/sqrt(2*n_max); %Standardavvikelse, medelvärde

%hold on
%plot(freq,10*log10(PSDx),'r')
%plot(freq,10*log10(PSDy),'m')
%hold off

%fBm

T=1e-2; %Tid mellan sampling
WVS_fGn=@(o,H)4*(sin(o*T/2)).^2.*abs(o.^(-(2*H+1)));
WVS_ejsin=@(o,H)T^2.*abs(o.^(-2*H+1));

o=freq;


o_bra=[(1:200) (227:250) (267:295) (317:380)]; %Tar ej med topparna
PSD_sum_o_bra=PSD_sum(o_bra);
freq_o_bra=freq(o_bra);

%H och normering för anpassad kurva
H_part=[0.2745 0.3287];
b_part=[17.43 20.76]/4;

b=b_part(fil);
H=H_part(fil);
H_ejsin=[0.2055 0.2886];
b_ejsin=[2.258 3.5];

hold on
%plot(o,10*log10(WVS_fGn(o,H)*b)) %Lite godtycklig normering
plot(o,WVS_fGn(o,H)*b)
plot(o,WVS_ejsin(o,H_ejsin(fil))*b_ejsin(fil)); %PSD, frekv
hold off



%Anpassat, cftool

%ed, utan peakar
% General model:
%      f(x) = a*(sin(1e-2*x/2)).^2.*x^(-(2*b+1))
% Coefficients (with 95% confidence bounds):
%        a =       1.607  (1.232, 1.982)
%        b =     0.00391  (-0.01846, 0.02628)
% 

%%%ed, innan peakar
% General model:
%      f(x) = a*(sin(1e-2*x/2)).^2.*x^(-(2*b+1))
% Coefficients (with 95% confidence bounds):
%        a =       17.43  (15.46, 19.41)
%        b =      0.2745  (0.2618, 0.2873)
% 

%%%ed, utan peakar ejsin
% General model:
%      f(x) = a*1e-4*x^(-2*b+1)
% Coefficients (with 95% confidence bounds):
%        a =       2.258  (1.957, 2.559)
%        b =      0.2055  (0.1929, 0.2182)


%%%%%%%%%%%%%%%%%%%%%

%log, alla utan peakar
% General model:
%      f(x) = a*(sin(1e-2/2*x)).^2.*x^(-(2*b+1))
% Coefficients (with 95% confidence bounds):
%        a =      0.7618  (0.586, 0.9375)
%        b =     0.00951  (-0.01261, 0.03162)
% 

%%%log, fram till peakar
% General model:
%      f(x) = a*(sin(1e-2*x/2)).^2.*x.^(-2*b-1)
% Coefficients (with 95% confidence bounds):
%        a =       20.76  (18.83, 22.7)
%        b =      0.3287  (0.3181, 0.3393)
% 

%%%log, utan peakar ejsin
%General model:
%     f(x) = a*1e-4*x^(-2*b+1)
% Coefficients (with 95% confidence bounds):
%        a =         3.5  (3.147, 3.853)
%        b =      0.2886  (0.2789, 0.2983)

%%

alpha=1.96; % 95% konfidensintervall
PSD_min=abs(PSD_sum-alpha*std_PSD);
PSD_max=abs(PSD_sum+alpha*std_PSD);

%H_max=[H_part(1)+0.10 H_part(2)+0.13]; %Övre gräns för H, ed först, log sen
%H_min=[H_part(1)-0.10 H_part(2)-0.13];
H_max=[0.2535 0.2789]; %ejsin
H_min=[0.1168 0.3657]; 

b_max=[0.40 0.5]; %Skalfaktor
b_min=[1.4 0.6];

H_maxsin=[0.3285 0.3451];
H_minsin=[0.1152 0.2663];

b_maxsin=[0.24 0.3];
b_minsin=[1.4 0.7];



hold on
plot(o,PSD_min,'r')
plot(o,PSD_max,'b')
plot(o,WVS_ejsin(o,H_max(fil)*b_max(fil)),'--b')
plot(o,WVS_ejsin(o,H_min(fil)*b_min(fil)),'--r')
plot(o,WVS_fGn(o,H_maxsin(fil)*b_maxsin(fil)),'b')
plot(o,WVS_fGn(o,H_minsin(fil)*b_minsin(fil)),'r')


hold off

%Anpassningar utan f=0, upp till f<147.1 rad/s

%%%ed, sin, halva
%max
% General model:
%      f(x) = a*(sin(1e-2*x/2)).^2.*x.^(-2*b-1)
% Coefficients (with 95% confidence bounds):
%        a =       39.84  (34.86, 44.83)
%        b =      0.3285  (0.3143, 0.3428)

%min
% Coefficients (with 95% confidence bounds):
%        a =       2.528  (2.144, 2.911)
%        b =      0.1152  (0.09864, 0.1317)


%%%%ed, ejsin, hela
%      f(x) = a*1e-4*x.^(-2*b+1)
% Coefficients (with 95% confidence bounds):
%        a =       4.916  (4.253, 5.578)
%        b =      0.2535  (0.2405, 0.2664)
%min
% Coefficients (with 95% confidence bounds):
%        a =      0.6013  (0.4958, 0.7068)
%        b =      0.1168  (0.1004, 0.1333)


%%%%%log, sin
%max
% General model:
%      f(x) = a*(sin(1e-2*x/2)).^2.*x.^(-2*b-1)
% Coefficients (with 95% confidence bounds):
%        a =       33.85  (30.25, 37.45)
%        b =      0.3451  (0.3329, 0.3573)
%min
% Coefficients (with 95% confidence bounds):
%        a =       6.984  (6.161, 7.807)
%        b =      0.2663  (0.2531, 0.2796)

%%%%%log, ejsin
%max
%      f(x) = a*1e-4*x.^(-2*b+1)
% Coefficients (with 95% confidence bounds):
%        a =       4.606  (4.01, 5.201)
%        b =      0.2789  (0.2665, 0.2914)

%min
% Coefficients (with 95% confidence bounds):
%        a =       3.914  (3.327, 4.501)
%        b =      0.3657  (0.3509, 0.3806)


%H_log=0.2384+0.09-0.13
%H_ed=0.2764+0.05-0.16

%%

% o_log=log10(o);
% PSD_max_log=log10(PSD_max);
% PSD_min_log=log10(PSD_min);
% PSD_log=log10(PSD_sum);

%ed
%max
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =      0.3678  (0.3593, 0.3763)
%        b =       -2.89  (-2.92, -2.86)
% 


%med
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =      0.3365  (0.3302, 0.3428)
%        b =      -3.161  (-3.183, -3.138)
% 




%min
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =      0.2273  (0.2173, 0.2374)
%        b =      -3.813  (-3.849, -3.778)
% 


%%%%%%
%log
%max
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =      0.3705  (0.3617, 0.3792)
%        b =      -3.226  (-3.258, -3.195)
% 


%med
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =      0.3363  (0.3299, 0.3428)
%        b =       -3.51  (-3.533, -3.487)
% 



%min
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =        0.22  (0.2087, 0.2313)
%        b =      -4.204  (-4.244, -4.164)
% 

b_max=[-2.89 -3.226]; %Skalfaktor
b_min=[-3.813 -4.204];
a_max=[0.3678 0.3705]; %Skalfaktor
a_min=[0.2273 0.22];

func=@(o,a,b)b*o.^(-2*a+1);


%plot(o,10*log10(WVS_fGn(o,H_max(fil))*b_max(fil)))
%plot(o,10*log10(WVS_fGn(o,H_min(fil))*b_min(fil)))


figure(4)
plot(o,PSD_sum)
hold on
plot(o,PSD_min)
plot(o,PSD_max)
plot(o,func(o,a_max(fil),b_max(fil)))
hold off
set(gca,'Xscale','log','Yscale','log')

%%
clc
plot_data=NaN(length(freq),6);
plot_data(:,1)=freq;
plot_data(:,2)=(PSD_sum);
plot_data(:,3)=(b*WVS_fGn(o,H));
plot_data(:,4)=(std_PSD);
plot_data(:,5)=

save('~/PSD_steg_ed.tsv', 'plot_data', '-ascii')

