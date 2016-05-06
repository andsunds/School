%PSD


clc;clf;clear all
load('filnamn.mat')
data_energydepleted =load(filnamn{1});
data_logphase=load(filnamn{2});

data=data_logphase;
C = separera(data);
n=length(C);%antal partiklar

Fs = 100;
n_max=915
t = C{1}(1:n_max,1);
x=diff(C{1}(1:n_max+1,2));

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;

plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

%%
%fBm

H=0.35;
T=0.01;
WVS_fGn=@(o)4*(sin(o*T/2)).^2.*abs(o.^(-(2*H+1)));

o=linspace(0,50,1000);
figure(5)
semilogy(o,WVS_fGn(o))