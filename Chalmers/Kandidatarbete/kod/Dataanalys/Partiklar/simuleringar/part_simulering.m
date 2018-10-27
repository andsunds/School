%%Brownian motion

%Kör autokorr.m, korrelation i N och T koordinater först för att få
%jämförelse mellan simulering och data

clear all;clc;%clf;
for j=1:10 %Bör köra fler för statistisk säkerhet men tar låg tid
hold on
n=304;% Antal partiklar (i energydepleted=304, logphase=193)
R=916; %Antal steg
X = randn(R,2,n);

XY = cumsum(X,1);
%figure(1)
%plot(XY(:,1),XY(:,2))

TN=zeros(1000,2);

%Antal "lags", dvs hur l???ngt i tiden vill man unders???ka korrelation
lagsauto =R-1;
lagscross = 900;
corrT = zeros(lagsauto+1,1);
corrN = zeros(lagsauto+1,1);


for i=1:n
    TN =  koordinatbyte( bsxfun(@minus, XY(:,:,i), mean(XY(:,:,i), 1)));
    T=TN(:,1);
    N=TN(:,2);
    corrT = corrT+autocorr(T,lagsauto); % Summera autokorrelation f???r alla partiklar
    corrN = corrN+autocorr(N,lagsauto); % Summera autokorrelation f???r alla partiklar
end
corrT=corrT/n;
corrN=corrN/n;

U = linspace(0,lagsauto,lagsauto+1); % Vektor med diskreta punkter, lags


figure(2)
subplot(1,2,1)
plot(U(1:end)*0.1,corrT(1:end),'*','MarkerSize',2)
axis([0 R*0.02 -0.2 1]) %För att samma plot som riktiga partiklarna
hold on
title('Korrelation T')
subplot(1,2,2)
plot(U(1:end)*0.1,corrN(1:end),'v','MarkerSize',2)
axis([0 R*0.02 -0.2 1])
hold on
title('Korrelation N')
end

%%

clc;clf;clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=2;
data =load(filnamn{fil});
C = separera(data);
n=length(C);%antal partiklar

gyr = zeros(2,2);
D = zeros(n,2,2);
V = zeros(n,2,2);
NT = cell(n,1);
for i=1:n
    X = C{i}(:,2)-mean(C{i}(:,2));
    Y = C{i}(:,3)-mean(C{i}(:,3));
    R = [C{i}(:,2),C{i}(:,3)];
    gyr(:,:) = [X'*X,X'*Y;
                Y'*X,Y'*Y]./length(C{i});
    [V(i,:,:),D(i,:,:)] = eig(gyr(:,:));
    Vort = [V(i,2,2),-V(i,1,2)];
    NT{i} = [R*V(i,:,2)' R*Vort'];  
end

%Antal "lags", dvs hur l???ngt i tiden vill man unders???ka korrelation
R=800;
lagsauto =R-1;
lagscross = 900;
corrT = zeros(lagsauto+1,1);
corrN = zeros(lagsauto+1,1);


for i=1:n 
    T=NT{i}(:,1);
    N=NT{i}(:,2);
    corrT = corrT+autocorr(T,lagsauto); % Summera autokorrelation f???r alla partiklar
    corrN = corrN+autocorr(N,lagsauto); % Summera autokorrelation f???r alla partiklar
end
corrT=corrT/n;
corrN=corrN/n;

U = linspace(0,lagsauto,lagsauto+1); % Vektor med diskreta punkter, lags


figure(2)
subplot(1,2,1)
plot(U(1:end),corrT(1:end),'*','MarkerSize',2)
hold on
title('Korrelation T')
subplot(1,2,2)
plot(U(1:end),corrN(1:end),'v','MarkerSize',2)
hold on
title('Korrelation N')






