%PSD

clc;
%clf;
clear all
load('filnamn.mat')
data_energydepleted =load(filnamn{1});
data_logphase=load(filnamn{2});
namn=[data_energydepleted data_logphase];

fil=2; %1=ed, 2=log

data=namn(fil);
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
plot(freq,10*log10(PSD_sum))
grid on
title('PSD för stegen')
xlabel('Frequency (rad/s)')
ylabel('Power/Frequency (dB/(rad/s))')

mean_PSD=mean(PSDx,2);
std_PSDx=std(PSDx,0,2);
std_PSDy=std(PSDy,0,2);
std_PSD=sqrt(std_PSDx.^2+std_PSDy.^2)/sqrt(2*n_max); %Standardavvikelse, medelvärde

%hold on
%plot(freq,10*log10(PSDx),'r')
%plot(freq,10*log10(PSDy),'m')
%hold off

%%
%fBm

T=1e-2; %Tid mellan sampling
WVS_fGn=@(o,H)4*(sin(o*T/2)).^2.*abs(o.^(-(2*H+1)));
WVS_ejsin=@(o,H)(T*o).^2.*abs(o.^(-(2*H+1)));

o=freq;

%H=0.22; %Startvärde
%index_norm=12;

%for i=1:20 %För att iterera fram H, pga normeringen
%norm_plot=((PSD_sum(index_norm)))/((WVS_fGn(freq(index_norm),H))); %För normering till data
%norm_plot=1;
%WVS_fGn_norm=@(o,H)4*(sin(o*T/2)).^2.*abs(o.^(-(2*H+1)))*norm_plot; %Normerade funktioner
%WVS_ejsin_norm=@(o,H)(T*o).^2.*abs(o.^(-(2*H+1)))*norm_plot;

o_bra=[(1:200) (227:250) (267:295) (317:380)]; %Tar ej med topparna
PSD_sum_o_bra=PSD_sum(o_bra);
freq_o_bra=freq(o_bra);
%o_bra=1:length(freq);
%diff_anp=@(H)abs((WVS_fGn_norm(o_bra,H)-PSD_sum(o_bra)'));

%H=lsqnonlin(diff_anp,0.1)

% hold on
% plot(o,10*log10(WVS_fGn_norm(o,H))) %Lite godtycklig normering
% plot(o,10*log10(WVS_ejsin_norm(o,H)),'-')
% xlabel('Frekvens (Hz)')
% ylabel('Effekt/Frekvens (dB/Hz)')
% hold off

%end


H_part=[0.28 0.24];
b_part=[4.6e0 1.5e0];

b=b_part(fil);
H=H_part(fil);

hold on
plot(o,10*log10(WVS_fGn(o,H)*b)) %Lite godtycklig normering
%plot(o,10*log10(WVS_ejsin(o,H)*b)); %PSD, frekv
hold off



%Anpassat, cftool

%ed, utan peakar
% General model:
%      f(x) = a*(sin(1e-2*x/2)).^2.*x^(-(2*b+1))
% Coefficients (with 95% confidence bounds):
%        a =       1.598  (1.231, 1.964)
%        b =    0.003314  (-0.01871, 0.02533)
% 
% Goodness of fit:
%   SSE: 0.00025
%   R-square: 0.7592
%   Adjusted R-square: 0.7586
%   RMSE: 0.0007955

%%%ed, innan peakar
% General model:
%      f(x) = a*(sin(1e-2*x/2)).^2.*x^(-(2*b+1))
% Coefficients (with 95% confidence bounds):
%        a =       17.62  (15.67, 19.57)
%        b =      0.2764  (0.264, 0.2889)
% 
% Goodness of fit:
%   SSE: 1.586e-05
%   R-square: 0.861
%   Adjusted R-square: 0.8603
%   RMSE: 0.0002755



%%%%%%%%%%

%log, alla utan peakar
% General model:
%      f(x) = a*(sin(1e-2/2*x)).^2.*x^(-(2*b+1))
% Coefficients (with 95% confidence bounds):
%        a =      0.7618  (0.586, 0.9375)
%        b =     0.00951  (-0.01261, 0.03162)
% 
% Goodness of fit:
%   SSE: 5.281e-05
%   R-square: 0.751
%   Adjusted R-square: 0.7503
%   RMSE: 0.0003638

%%%log, fram till peakar
% General model:
%      f(x) = a*(sin(1e-2*x/2)).^2.*x.^(-2*b-1)
% Coefficients (with 95% confidence bounds):
%        a =       20.17  (18.26, 22.08)
%        b =      0.3259  (0.3151, 0.3367)
% 
% Goodness of fit:
%   SSE: 8.394e-06
%   R-square: 0.7918
%   Adjusted R-square: 0.7908
%   RMSE: 0.000199

%%

alpha=1.96; % 95% konfidensintervall
PSD_min=abs(PSD_sum-alpha*std_PSD);
PSD_max=abs(PSD_sum+alpha*std_PSD);

H_max=[H_part(1)+0.10 H_part(2)+0.13]; %Övre gräns för H, ed först, log sen
H_min=[H_part(1)-0.10 H_part(2)-0.13];
b_max=[1.5e1 7e0]; %Skalfaktor
b_min=[11e-1 3e-1];



hold on
plot(o,10*log10(PSD_min))
plot(o,10*log10(PSD_max))
plot(o,10*log10(WVS_fGn(o,H_max(fil))*b_max(fil)))
plot(o,10*log10(WVS_fGn(o,H_min(fil))*b_min(fil)))
hold off

%Anpassningar utan f=0, upp till f<147.1 rad/s

%%%%ed
%max
% Coefficients (with 95% confidence bounds):
%        a =       39.63  (34.7, 44.56)
%        b =      0.3285  (0.3143, 0.3426)
% 
% Goodness of fit:
%   SSE: 5.417e-05
%   R-square: 0.704
%   Adjusted R-square: 0.7026
%   RMSE: 0.0005055

%min
% Coefficients (with 95% confidence bounds):
%        a =       2.578  (2.189, 2.967)
%        b =      0.1175  (0.101, 0.1339)
% 
% Goodness of fit:
%   SSE: 6.109e-06
%   R-square: 0.9338
%   Adjusted R-square: 0.9335
%   RMSE: 0.0001698

%%%%%log
%max
% Coefficients (with 95% confidence bounds):
%        a =       17.71  (15.35, 20.07)
%        b =        0.33  (0.3148, 0.3452)
% 
% Goodness of fit:
%   SSE: 1.218e-05
%   R-square: 0.6782
%   Adjusted R-square: 0.6767
%   RMSE: 0.0002397

%min
% Coefficients (with 95% confidence bounds):
%        a =       1.019  (0.8581, 1.179)
%        b =      0.1066  (0.08941, 0.1238)
% 
% Goodness of fit:
%   SSE: 1.226e-06
%   R-square: 0.9326
%   Adjusted R-square: 0.9322
%   RMSE: 7.604e-05

%H_log=0.2384+0.09-0.13
%H_ed=0.2764+0.05-0.16

%%

o_log=log10(o);
PSD_max_log=log10(PSD_max);
PSD_min_log=log10(PSD_min);
PSD_log=log10(PSD_sum);

%ed
%max
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =      0.3678  (0.3593, 0.3763)
%        b =       -2.89  (-2.92, -2.86)
% 
% Goodness of fit:
%   SSE: 0.5866
%   R-square: 0.8156
%   Adjusted R-square: 0.8148
%   RMSE: 0.05235

%med
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =      0.3365  (0.3302, 0.3428)
%        b =      -3.161  (-3.183, -3.138)
% 
% Goodness of fit:
%   SSE: 0.2995
%   R-square: 0.9271
%   Adjusted R-square: 0.9268
%   RMSE: 0.03813



%min
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =      0.2273  (0.2173, 0.2374)
%        b =      -3.813  (-3.849, -3.778)
% 
% Goodness of fit:
%   SSE: 0.7612
%   R-square: 0.933
%   Adjusted R-square: 0.9327
%   RMSE: 0.06079

%%%%%%
%log
%max
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =      0.3705  (0.3617, 0.3792)
%        b =      -3.226  (-3.258, -3.195)
% 
% Goodness of fit:
%   SSE: 0.5788
%   R-square: 0.8052
%   Adjusted R-square: 0.8043
%   RMSE: 0.05301

%med
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =      0.3363  (0.3299, 0.3428)
%        b =       -3.51  (-3.533, -3.487)
% 
% Goodness of fit:
%   SSE: 0.3154
%   R-square: 0.9237
%   Adjusted R-square: 0.9234
%   RMSE: 0.03913


%min
% General model:
%      f(x) = (-2*a+1)*x+b
% Coefficients (with 95% confidence bounds):
%        a =        0.22  (0.2087, 0.2313)
%        b =      -4.204  (-4.244, -4.164)
% 
% Goodness of fit:
%   SSE: 0.9743
%   R-square: 0.9203
%   Adjusted R-square: 0.9199
%   RMSE: 0.0686

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

