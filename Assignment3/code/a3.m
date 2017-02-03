%% 
% U^{n+1}_j=(1-beta)*U_{j}+(beta-alpha)/2*U_{j+1}+(alpha-beta)/2*U_{j-1}
clc;clear;clf

J=1001;
N=100;
a=1;
Dx=2/(J-1); %we have to have (J-1) here, to get J values in x.
x=(-1:Dx:1).';
alpha=.99;
Dt=alpha*Dx/a;


%different beta for different schemes.
% scheme |   LW    | LF |   upwind   | centered |
% beta   | alpha^2 | 1  | abs(alpha) |    0     |
beta=alpha^2;


T=spTranferMatrix(J,alpha,beta);

k=10;
IC=sin(2*pi*k*x);

U=zeros(J,N);
U(:,1)=IC;

for n=1:N-1;
    U(:,n+1)=T*U(:,n);
end

%
y=U(:,1);
p=plot(x,y);
p.YDataSource='y';

for i=2:N;
    pause(.3)
    y=U(:,i);
    refreshdata
    title(sprintf('n = %d',i))
end