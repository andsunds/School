%CTRW

%Simulera CTRW
n=1;% Antal partiklar
R=900000;%Anpassas efter hur alpha v?ljs

%Generera v?ntetider, lite hokus pokus och sopande under mattan
F_i=@(U,eps,alpha,A)(eps^(-alpha)-alpha*U/A).^(-1/alpha); %Invers CDF, CDF=tau^-(1+alpha)
alpha=0.70; %Anpassningsbar parameter
eps=1e-9; %Minsta m?jliga v?ntetid
A=alpha*eps^alpha; %Normeringskonstant
U=rand(R,n);
tau=F_i(U,eps,alpha,A); %V?ntetider
t=cumsum(tau,1); %tid

X = randn(R,2,n)*1e-9;%Normalf?rdelade steg

XY=cumsum(X,1);
plot(t,XY)




