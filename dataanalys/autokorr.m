% Autokorrelation samt "cross-correlation"
clear all;clc;clf;

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=1;
data = separera(load(filnamn{fil}));

%Antal "lags", dvs hur långt i tiden vill man undersöka korrelation
lagsauto =900;
lagscross = 900;
corr = zeros(lagsauto+1,1);
cross = zeros(2*lagscross+1,1); 

%% En partikel
i =100;%Välj partikel
XY = data{i}(:,2:3);
avst1=sqrt(sum(XY.^2,2)); 

figure(1)
autocorr(avst1) %Korrelation av storlek på steg

%% Alla partiklar 

for i=1:length(data)
    if data{i}(:,4)<20 % Gräns för storlek på partiklar
    XY = data{i}(:,2:3);
    avst=sqrt(sum(XY.^2,2)); 
    corr = corr+autocorr(avst,lagsauto); % Summera autokorrelation för alla partiklar
    cross = cross+crosscorr(diff(data{i}(:,2)),diff(data{i}(:,3)),lagscross);
    end
end
corr = corr/length(data);%Normera med antalet partiklar
cross = cross/length(data);

U = linspace(0,lagsauto,lagsauto+1); % Vektor med diskreta punkter, lags
V = linspace(-lagscross+1,lagscross+1,2*lagscross+1); % Vektor med diskreta punkter, lags

figure(2)
subplot(1,2,1)
plot(U(2:end),corr(2:end),'*')
title('Autokorrelation')
subplot(1,2,2)
plot(V,cross,'v')
title('Crosskorrelation')


%% Korrelation i normal och tangential koordinater

%Antal "lags", dvs hur långt i tiden vill man undersöka korrelation
lagsauto =900;
lagscross = 900;
corrT = zeros(lagsauto+1,1);
corrN = zeros(lagsauto+1,1);
for i=1:length(data)
    if data{i}(:,4)<20 % Gräns för storlek på partiklar
    
    TN=koordinatbyte( bsxfun(@minus, data{i}(:,2:3), mean(data{i}(:,2:3), 1)) );
    T=TN(:,1);
    N=TN(:,2);
    corrT = corrT+autocorr(T,lagsauto); % Summera autokorrelation för alla partiklar
    corrN = corrN+autocorr(N,lagsauto); % Summera autokorrelation för alla partiklar
    
    end
end
corrT=corrT/length(data);
corrN=corrN/length(data)

U = linspace(0,lagsauto,lagsauto+1); % Vektor med diskreta punkter, lags

figure(2)
subplot(1,2,1)
plot(U(2:end),corrT(2:end),'*')
title('Korrelation T')
subplot(1,2,2)
plot(U,corrN,'v')
title('Korrelation N')
