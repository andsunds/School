% Autokorrelation samt "cross-correlation"
clear all;clc;clf;

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=1;
data = separera(load(filnamn{fil}));

%Antal "lags", dvs hur l???ngt i tiden vill man unders???ka korrelation
lagsauto =900;
lagscross = 900;
corr = zeros(lagsauto+1,1);
cross = zeros(2*lagscross+1,1); 

%% En partikel
clf

i =100;%Choose particle
XY = data{i}(:,2:3);
avst1=sqrt(sum(XY.^2,2)); 

figure(1)
autocorr(avst1) %Correlation of step length

%% Alla partiklar 
clf

for i=1:length(data)
    if data{i}(:,4)<20 % Upper limit for particle size
    XY = data{i}(:,2:3);
    avst=sqrt(sum(XY.^2,2)); 
    corr = corr+autocorr(avst,lagsauto); % Sum of autocorrelation for all particles
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
axis([V(1) V(end) -0.08 0.08])
title('Korskorrelation')


%% Korrelation i normal och tangential koordinater
clf

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

for j=1:2
    
fil=j;

%Antal "lags", dvs hur l???ngt i tiden vill man unders???ka korrelation
lagsauto =900;
lagscross = 900;
corrT = zeros(lagsauto+1,1);
corrN = zeros(lagsauto+1,1);

data = separera(load(filnamn{fil}));
for i=1:length(data)
    if data{i}(1,4)<40 % Gr???ns f???r storlek p??? partiklar
    
    TN=koordinatbyte( bsxfun(@minus, data{i}(:,2:3), mean(data{i}(:,2:3), 1)) );
    T=TN(:,1);
    N=TN(:,2);
    corrT = corrT+autocorr(T,lagsauto); % Summera autokorrelation f???r alla partiklar
    corrN = corrN+autocorr(N,lagsauto); % Summera autokorrelation f???r alla partiklar
    
    end
end
corrT=corrT/length(data);
corrN=corrN/length(data);

U = linspace(0,lagsauto,lagsauto+1); % Vektor med diskreta punkter, lags


figure(2)
subplot(1,2,1)
U_s=U*0.1; %To plot in seconds if 10fps
plot(U_s,corrT,'*','MarkerSize',3)
set(gca,'fontsize',16)
hold on
title('Korrelation t-koordinat','fontsize',18)
subplot(1,2,2)
plot(U_s,corrN,'*','MarkerSize',3)
set(gca,'fontsize',16)
hold on
title('Korrelation n-koordinat','fontsize',18)


end

subplot(1,2,1)
legend('Dvala','Aktiva','Location','Best') %Placerar en legend p?? l??mplig plats
axis([min(U_s) 20 -0.1 1])
xlabel('dt (s)','fontsize',16)
subplot(1,2,2)
legend('Dvala','Aktiva','Location','Best')
axis([min(U_s) 20 -0.1 1])
xlabel('dt (s)','fontsize',16)

%%
%Korskorrelation i N och T-variablerna
%Hur de båda signalerna liknar varandra

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

data = separera(load(filnamn{fil}));

lagscross = 900;
cross = zeros(2*lagscross+1,1); 

for i=1:length(data)
    %if data{i}(:,4)<40 % Upper limit for particle size
    TN=koordinatbyte( bsxfun(@minus, data{i}(:,2:3), mean(data{i}(:,2:3), 1)) );
    T=TN(:,1);
    N=TN(:,2);
    %cross = cross+crosscorr(diff(data{i}(:,2)),diff(data{i}(:,3)),lagscross);
    cross = cross+crosscorr(diff(T),diff(N),lagscross);
    %end
end

cross = cross/length(data);
V = linspace(-lagscross+1,lagscross+1,2*lagscross+1); % Vektor med diskreta punkter, lags

figure(3)
plot(V,cross,'v')
axis([V(1) V(end) -0.08 0.08])
title('Korskorrelation')
