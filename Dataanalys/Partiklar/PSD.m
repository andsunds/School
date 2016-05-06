%PSD


clc;clf;clear all
load('filnamn.mat')
data_energydepleted =load(filnamn{1});
data_logphase=load(filnamn{2});

data=data_energydepleted;
C = separera(data);
n=length(C);%antal partiklar

Fs = 100;
n_max=915
t = C{1}(1:n_max,1);
PSD=0;

for i=1:n

x=diff(C{i}(1:n_max+1,2));
%x=

N = length(x);
xdft = fft(x); %Fouriertransformera
xdft = xdft(1:N/2+1); %Ta bara halva spektrumet
psdx = (1/(Fs*N)) * abs(xdft).^2; %PSD
psdx(2:end-1) = 2*psdx(2:end-1); %För att bevara totala effekten, vi tog ju bort halva ovan
%0 och Nyquist frek förekommer bara 1 gång därav ej *2

PSD=PSD+psdx;
end

PSD=PSD/n_max;
freq = 0:Fs/length(x):Fs/2;
%freq = 0:(2*pi)/N:pi; %Vinkelfrekv

figure(5)
plot(freq,10*log10(PSD))
grid on
title('PSD för stegen')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

%%
%fBm

H=0.28;
T=1e-2;
WVS_fGn=@(o)4*(sin(o*T/2)).^2.*abs(o.^(-(2*H+1)));
WVS_gauss=@(o)abs(o.^(-(2*0.5+1)));

%o=linspace(0,50,1000);
o=freq;

hold on
plot(o,10*log10(WVS_fGn(o)*3.5e-15)) %Lite godtycklig normering
%plot(o,10*log10(WVS_gauss(o)),'-')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
hold off