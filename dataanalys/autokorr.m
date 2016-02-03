% Autokorrelation samt "cross-correlation"
clear all;clc;clf;

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=1;
data = separera(load(filnamn{fil}));

%Antal "lags", dvs hur l�ngt i tiden vill man unders�ka korrelation
lagsauto =900;
lagscross = 900;
corr = zeros(lagsauto+1,1);
cross = zeros(2*lagscross+1,1); 

%% En partikel
clf

i =100;%V�lj partikel
XY = data{i}(:,2:3);
avst1=sqrt(sum(XY.^2,2)); 

figure(1)
autocorr(avst1) %Korrelation av storlek p� steg

%% Alla partiklar 

for i=1:length(data)
    if data{i}(:,4)<20 % Gr�ns f�r storlek p� partiklar
    XY = data{i}(:,2:3);
    avst=sqrt(sum(XY.^2,2)); 
    corr = corr+autocorr(avst,lagsauto); % Summera autokorrelation f�r alla partiklar
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
clf

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

for j=1:2
    
fil=j;

%Antal "lags", dvs hur l�ngt i tiden vill man unders�ka korrelation
lagsauto =900;
lagscross = 900;
corrT = zeros(lagsauto+1,1);
corrN = zeros(lagsauto+1,1);

data = separera(load(filnamn{fil}));
for i=1:length(data)
    if data{i}(1,4)<22 % Gr�ns f�r storlek p� partiklar
    
    TN=koordinatbyte( bsxfun(@minus, data{i}(:,2:3), mean(data{i}(:,2:3), 1)) );
    T=TN(:,1);
    N=TN(:,2);
    corrT = corrT+autocorr(T,lagsauto); % Summera autokorrelation f�r alla partiklar
    corrN = corrN+autocorr(N,lagsauto); % Summera autokorrelation f�r alla partiklar
    
    end
end
corrT=corrT/length(data);
corrN=corrN/length(data);

U = linspace(0,lagsauto,lagsauto+1); % Vektor med diskreta punkter, lags


figure(2)
subplot(1,2,1)
plot(U(2:end),corrT(2:end),'*','MarkerSize',2)
hold on
title('Korrelation T')
subplot(1,2,2)
plot(U,corrN,'v','MarkerSize',2)
hold on
title('Korrelation N')


end

subplot(1,2,1)
legend('Dvala','Aktiva','Location','Best') %Placerar en legend på lämplig plats
subplot(1,2,2)
legend('Dvala','Aktiva','Location','Best')