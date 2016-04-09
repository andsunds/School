%CTRW

%F?rs?k till att simulera CTRW
n=1;% Antal partiklar
R=9000000;%Anpassas efter hur alpha v?ljs

%Generera v?ntetider, lite hokus pokus och sopande under mattan
F_i=@(U,eps,alpha,A)(eps^(-alpha)-alpha*U/A).^(-1/alpha); %Invers CDF, CDF=tau^-(1+alpha)
alpha=0.85; %Anpassningsbar parameter
eps=1e-9; %Minsta m?jliga v?ntetid
A=alpha*eps^alpha; %Normeringskonstant
U=rand(R,n);
tau=F_i(U,eps,alpha,A); %V?ntetider
t=cumsum(tau,1); %tid

X = randn(R,2,n);%Normalf?rdelade steg

XY=cumsum(X,1); %x,y separerade
R=sqrt(XY(:,1).^2+XY(:,2).^2);
subplot(2,1,1)
plot(t,XY)
xlabel('tid (s)')
ylabel('F?rflyttning')
title('XY')
subplot(2,1,2)
plot(t,R)
xlabel('tid (s)')
ylabel('F?rflyttning')
title('R')



