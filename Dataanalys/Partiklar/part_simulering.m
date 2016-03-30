%% 

clear all;clc;
hold on
n=1000;% Antal partiklar
R=3000;
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
plot(U(1:end),corrT(1:end),'*','MarkerSize',2)
hold on
title('Korrelation T')
subplot(1,2,2)
plot(U(1:end),corrN(1:end),'v','MarkerSize',2)
hold on
title('Korrelation N')




