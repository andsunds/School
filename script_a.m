%%
clc; clf;clear all
p='Data/';addpath(p);

files_t=dir([p,'Pb','*_t.lvm']);
files_v=dir([p,'Pb','*_v.lvm']);

l=length(files_t);

% Antal övertoner+grundton
m = 5;
T=zeros(l,m);
Pb_order = [2 3 4 6 5]-1;
Pb_x = [0,1,2,3,4].';
betaC  = zeros(l,m); betaV  = zeros(l,m);
gammaC = zeros(l,m); gammaV = zeros(l,m);

for i=1:l
    name_t=files_t(i).name;
    name_v=files_v(i).name;
    time = load(name_t);
    % Konvertera till sekunder
    time = (time-min(min(time)))/1000;
    volt = load(name_v);
    % Konvertera till temperatur
    cels = V2C(volt,Pb_order);
    % Läs in excitationens periodtid
    T(i,1) = str2double(name_t(regexp(name_t, '\d')));
    sample_count = size(time,1);
    
    % Allokera minne för en variabel som den återskapade signalen placeras
    % i efter att frekvenskomponenterna beräknas med heterodynmetoden. Den
    % initialiseras med DC-värdet (pga transienter förväntas en offset).
    % Alternativt skulle man kunna ha ett löpande medelvärde som är två
    % perioder långt.
    %num_points = timT(i,1)/(time(end)-time(1))
    cels_hat = zeros(size(time,1),5);
    for j=1:5
        cels_hat(:,j) = smooth(time(:,j),cels(:,j), 60*T(i,1),'moving');
    end
    
    % Högpassfiltrera?
    %cels = highpass(time,cels,0.1/T(j));
    
    subplot(2,2,[1,2])
    plot(repmat(time(:,1)/60,1,5),cels)
    title(sprintf('Periodtid: %dmin',T(i)))
    xlabel('Tid [min]')
    ylabel('Temperatur [C]')
    
    for j=1:m
        Nper = 5;
         % Se till så att ett helt antal perioder av den exciterande
         % vågformen kommer med
        T(i,j) = T(i,1)/j;
        Nper = Nper*j;
        % Plocka ut amplitud och fas för denna frekvenskomponent
        [log_rC, phiC] = Heterodyn(time,cels,T(i,j)*60, Nper, Pb_order);
        [log_rV, phiV] = Heterodyn(time,volt,T(i,j)*60, Nper, Pb_order);
        
        %%%%%%Roetzels metod%%%%%%
        %(phi1(:)-phi1(1))
        %k=-(phi1(1:N)-phi1(1))\(log_r1(1:N)-log_r1(1))
        %k=mean((log_r1(2:4)-log_r1(1))./(phi1(2:4)-phi1(1)))
        %tau=((1-k.^2)./(2*k))*T(i,j)/2/pi
        
        % Antal punkter i den linjära regressionen
        N=5;
        
        % Linjär regression för att beräkna faskonstant och dämpning
        XC = lscov([ones(N,1),Pb_x(1:N)],[log_rC(1:N),phiC(1:N)],exp(log_rC));
        log_V0 = XC(1,1);  gamma = -XC(2,1);
        phi0   = XC(1,2);  beta  = -XC(2,2);
        % Spara undan koefficienter
        gammaC(i,j) = -XC(2,1);
        betaC(i,j) = -XC(2,2);
        
        XV = [ones(N,1),Pb_x(1:N)]\[log_rV(1:N),phiV(1:N)];
        gammaV(i,j) = -XV(2,1);
        betaV(i,j) = -XV(2,2);
        
        % Lägg till denna frekvenskomponent i den rekonstruerade signalen
        w = 2*pi/T(i,j);
        t = repmat(time(:,1)/60,1,5);
        x = repmat(Pb_x(Pb_order)',size(time,1),1);
        % Använd regressionskoefficienter
        cels_hat = cels_hat + cos(phi0+w*t-beta*x).*exp(log_V0-gamma*x);
        phiC = repmat(phiC', sample_count,1);
        log_rC = repmat(log_rC', sample_count,1);
        % Använd frekvensvärdena direkt
        %cels_hat = cels_hat + cos(phiC+w*t).*exp(log_rC);
        
        % Plotta mätt signal tillsammans med rekonstruerad signal
        subplot(2,2,[1,2])
        plot(repmat(time(:,1)/60,1,5),cels)
        legend(reshape(sprintf('%d cm',Pb_x(Pb_order)),4,5)')
        title(sprintf('Periodtid: %dmin, #Frekvenskomponenter: %d',T(i),j))
        xlabel('Tid [min]')
        ylabel('Temperatur [C]')
        hold on
        plot(t(:,1:2),cels_hat(:,Pb_order(1:2)),'--')
        hold off
    
        % Plotta regressionen
        subplot(2,2,3)
        plot(Pb_x,log_rC,'o')
        hold on
        plot(Pb_x, log_V0 - Pb_x*gammaC(i,j),'-')
        hold off
        subplot(2,2,4)
        plot(Pb_x,phiC,'v')
        hold on
        plot(Pb_x, phi0 - Pb_x*betaC(i,j),'-')
        hold off
        
        pause(0.1)
       
    end
end

%% Plotta dispersionsrelationen
% Lägg all data i endimensionella vektorer
w = 2*pi./reshape(T,m*l,1); % Vinkelfrekvens
betaC = reshape(betaC,m*l,1);
gammaC = reshape(gammaC,m*l,1);
clf

plot(w,gammaC,'o',w,betaC,'d')
legend('\beta, Realdel','\gamma, Imaginärdel')
xlabel('\omega/[rad/min]')
ylabel('k/[cm^{-1}]')
grid on
hold on
% Plotta vilka mätpunkter som hör ihop (H för Harmonic)
H = 5;
for i=(1+(H-1)*l):(H*size(T,1))
    plot([w(i),w(i)],[gammaC(i),betaC(i)])
end
hold off

%%
clf
w = 2*pi*T.^(-1);

a=3.4e-1;
tau=.01;
x=linspace(0,4);

BETA=@(tau, a) sqrt(w/(2*a)).*sqrt(sqrt((tau*w).^2+1)+tau*w);

GAMMA=@(tau,a) sqrt(w/(2*a)).*sqrt(sqrt((tau*w).^2+1)-tau*w);

X=fminsearch(@(X) sum( (BETA(X(1), X(2)) - betaC).^2 ), [1, 1])
Y=fminsearch(@(X) sum( (GAMMA(X(1), X(2)) - gammaC).^2 ), [1, 1])

tau=X(1);a=X(2);
y1=sqrt(x/(2*a)).*sqrt(sqrt((tau*x).^2+1)+tau*x);
tau=Y(1);a=Y(2);
y2=sqrt(x/(2*a)).*sqrt(sqrt((tau*x).^2+1)-tau*x);


subplot(1,2,1)
plot(w,betaC,'bo'),hold on
plot(w,betaV,'rx')
plot(x,y1)
grid on
hold off

title('beta')
legend('temperatur', 'spänning', 'location', 'Best')





subplot(1,2,2)
plot(w,gammaC,'bv'), hold on
plot(w,gammaV,'r^')
plot(x,y2)
grid on
hold off
%axis([0,max(tau)*1.1,0,max(gamma1)*1.1])
title('gamma')
legend('temperatur', 'spänning', 'location', 'best')

%%
pause()

clc;clf

plot(w,betaC./gammaC,'bo'),hold on
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

