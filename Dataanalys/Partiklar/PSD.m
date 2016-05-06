%PSD

clc;clf;clear all
load('filnamn.mat')
data_energydepleted =load(filnamn{1});
data_logphase=load(filnamn{2});

data=data_energydepleted;
C = separera(data);
n=length(C);%antal partiklar

Fs = 100;
n_max=914
t = C{1}(1:n_max,1);
PSDx=zeros(458,n);
PSDy=zeros(458,n);
fil=1;

load('filnamn.mat')
load(kompl,...
      'intensitet', 'medelsteg', 'sigma_brus', '-mat')

I=intensitet{fil};

%Detta är just nu samma normering som för rörligheten...
lambda=medelsteg{fil}-0e-9; %Räknar ej med bakgrundsbrus
%beräkna koefficienter för lutning i loglog
koef=[ones(size(I)), log(I)]\log(lambda);
koef(1)=exp(koef(1));%konverterar till potenssamband

% plot(I,medelsteg{fil}-1e-9,'o')
% x=logspace(-1,1.3);
% y=koef(1)*x.^(koef(2));
% hold on
% plot(x,y)
% set(gca,'xscale','log','yscale','log')


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
psdx = (1/(Fs*N)) * abs(xdft).^2; %PSD
psdx(2:end-1) = 2*psdx(2:end-1); %För att bevara totala effekten, vi tog ju bort halva ovan
%0 och Nyquist frek förekommer bara 1 gång därav ej *2

ydft = fft(dy); %Fouriertransformera
ydft = ydft(1:N/2+1); %Ta bara halva spektrumet
psdy = (1/(Fs*N)) * abs(ydft).^2; %PSD
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

PSD_sum=(sum(PSDx,2)+sum(PSDy,2))/2;


freq = 0:Fs/length(dx):Fs/2;
%freq = 0:(2*pi)/N:pi; %Vinkelfrekv

figure(5)
plot(freq,10*log10(PSD_sum))
grid on
title('PSD för stegen')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

mean_PSD=mean(PSDx,2);
std_PSD=std(PSDx');

hold on
plot(freq,10*log10(PSD_sum)+std_PSD','-k')
plot(freq,10*log10(PSD_sum)-std_PSD','-k')
hold off
%%
%fBm

H=0.28;
T=1e-2;
WVS_fGn=@(o)4*(sin(o*T/2)).^2.*abs(o.^(-(2*H+1)));
WVS_ejsin=@(o)(T*o).^2.*abs(o.^(-(2*H+1)));

%o=linspace(0,50,1000);
o=freq;

hold on
plot(o,10*log10(WVS_fGn(o)*8e0)) %Lite godtycklig normering
plot(o,10*log10(WVS_ejsin(o)*8e0),'-')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
hold off