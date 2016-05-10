%PSD

clc;clf;clear all
load('filnamn.mat')
data_energydepleted =load(filnamn{1});
data_logphase=load(filnamn{2});

data=data_energydepleted;
C = separera(data);
n=length(C);%antal partiklar

Fs = 100;
n_max=914;
t = C{1}(1:n_max,1);
PSDx=zeros(458,n);
PSDy=zeros(458,n);

fil=2; %1=ed, 2=log

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

PSDx=PSDx/n_max;
PSDy=PSDy/n_max;

PSD_sum=(sum(PSDx,2)+sum(PSDy,2))/2; %Totala PSD:n


freq = (0:Fs/length(dx):Fs/2)/(2*pi); %Vinkelfrekvens läng x-axeln
%freq = 0:(2*pi)/N:pi; %Vinkelfrekv, 0-pi

figure(5)
plot(freq,10*log10(PSD_sum))
grid on
title('PSD för stegen')
xlabel('Frequency (rad/s)')
ylabel('Power/Frequency (dB/(rad/s))')

mean_PSD=mean(PSDx,2);
std_PSDx=std(PSDx');
std_PSDy=std(PSDy');
std_PSD=sqrt(std_PSDx.^2+std_PSDy.^2); %Standardavvikelse

%hold on
%plot(freq,10*log10(PSDx),'r')
%plot(freq,10*log10(PSDy),'m')
%hold off

%fBm

T=1e-2; %Tid mellan sampling
WVS_fGn=@(o,H)4*(sin(o*T/2)).^2.*abs(o.^(-(2*H+1)));
WVS_ejsin=@(o,H)(T*o).^2.*abs(o.^(-(2*H+1)));

o=freq;

H=0.22; %Startvärde
index_norm=12;

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

H_ed=0.2336;
H_log=0.2356;
b_ed=2.1e1;
b_log=0.9e1;
b=b_ed;
H=H_ed;
hold on
plot(o,10*log10(b*WVS_fGn(o,H))) %Lite godtycklig normering
%plot(o,10*log10(WVS_ejsin_norm(o,H))); %PSD, frekv
hold off



%Anpassat, cftool

%ed
% General model:
%      f(x) = b*(x).^(-2*a+1)
% Coefficients (with 95% confidence bounds):
%        a =      0.2336  (0.2226, 0.2446)
%        b =    0.002068  (0.002002, 0.002135)
% 
% Goodness of fit:
%   SSE: 5.053e-05
%   R-square: 0.9192
%   Adjusted R-square: 0.9189
%   RMSE: 0.0004031


%log
% General model:
%      f(x) = b*(x).^(-2*a+1)
% Coefficients (with 95% confidence bounds):
%        a =      0.2356  (0.2249, 0.2463)
%        b =   0.0009274  (0.0008984, 0.0009565)
% 
% Goodness of fit:
%   SSE: 9.555e-06
%   R-square: 0.9214
%   Adjusted R-square: 0.9212
%   RMSE: 0.0001756

%%
clc
plot_data=NaN(length(freq),4);
plot_data(:,1)=freq;
plot_data(:,2)=10*log10(PSD_sum);
plot_data(:,3)=10*log10(WVS_fGn_norm(o,H));
plot_data(:,4)=10*log10(std_PSD);

save('~/PSD_steg_log.tsv', 'plot_data', '-ascii')

