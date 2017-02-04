%% Playing around
% U^{n+1}_j=(1-beta)*U_{j}+(beta-alpha)/2*U_{j+1}+(alpha-beta)/2*U_{j-1}
clc;clear;clf

J=1001;
a=1;
Dx=2/(J-1); %we have to have (J-1) here, to get J values in x.
x=(-1:Dx:1).';
alpha=.2;
Dt=alpha*Dx/a;
T=0.1;
N=floor(T/Dt);


%different beta for different schemes.
% scheme |   LW    | LF |   upwind   | centered |
% beta   | alpha^2 | 1  | abs(alpha) |    0     |
beta=alpha^2;


T=spTranferMatrix(J,alpha,beta);

%k1=50;
%k2=2;
%IC=sin(2*pi*k1*x)+sin(2*pi*k2*x);
sigma=.001;
IC=exp(-x.^2/(2*sigma));


U=zeros(J,N);
U(:,1)=IC;

for n=1:N-1;
    U(:,n+1)=T*U(:,n);
end

%% Movie

y=U(:,1);
p=plot(x,y);
axis([-1,1, -1-min(U(:,1)),max(U(:,1))])
p.YDataSource='y';
%set(gca, 'yscale','log')

%I=find(max(U,[],1)==U(floor(J/2),:));

for i=1:1:N;
    pause(.01)
    y=U(:,i);
    refreshdata
    title(sprintf('t = %0.0f * 10^{-3}',i*Dt*1e3), 'fontsize', 20)
end


%% real work 1a-c
% U^{n+1}_j=(1-beta)*U_{j}+(beta-alpha)/2*U_{j+1}+(alpha-beta)/2*U_{j-1}
clc;clear;clf

lines={'-k',':k','-.k','--k'};

%JJ=[501,141,101]
%J=JJ(j);

J=101;
a=1;
Dx=2/(J-1); %we have to have (J-1) here, to get J values in x.
x=(-1:Dx:1).';

%ALPHA=[0.1,0.5,0.9];

%for j=1:3
alpha=.6;
%alpha=ALPHA(j);

Dt=alpha*Dx/a;
t=10;
N=floor(t/Dt);

%different beta for different schemes.
% scheme |   LW    | LF |   upwind   | centered |
% beta   | alpha^2 | 1  | abs(alpha) |    0     |


%BETA=[alpha^2, 1, abs(alpha), 0];

%for j=1:4
%beta=BETA(j);
beta=alpha^2;

T=spTranferMatrix(J,alpha,beta);

%k1=3;
%IC=sin(2*pi*k1*x);
sigma=0.01;
IC=exp(-x.^2/(2*sigma));


tic
TN=T^N;
%[V,D]=eigs(T,J);
%TN=V*(D.^N)/V;
toc

plot(x,TN*IC,lines{j},'linewidth',3)
hold on
%end
%%
%L=legend('LW',...
%         'Upwind','Centered');
L=legend('$\Delta{t}=2\times10^{-3}$',...
         '$\Delta{t}=10\times10^{-3}$','$\Delta{t}=18\times10^{-3}$');
set(L,'fontsize',30,'interpreter','latex')

set(gca,'fontsize',30, 'ylim',[-1,1])




%% real work 1d
% U^{n+1}_j=(1-beta)*U_{j}+(beta-alpha)/2*U_{j+1}+(alpha-beta)/2*U_{j-1}
clc;clear;clf


J=101;
a=1;
Dx=2/(J-1); %we have to have (J-1) here, to get J values in x.
x=(-1:Dx:1).';


alpha=.6;


Dt=alpha*Dx/a;
t=10.2;
N=floor(t/Dt);

%different beta for different schemes.
% scheme |   LW    | LF |   upwind   | centered |
% beta   | alpha^2 | 1  | abs(alpha) |    0     |


beta=alpha^2;

T=spTranferMatrix(J,alpha,beta);

sigma=0.005;
IC=exp(-x.^2/(2*sigma));


tic
TN=T^N;
toc

plot(x,TN*IC,'-k','linewidth',3)
hold on
plot(x,IC,':k','linewidth',3)

L=legend('$U(x, t_{0}=0)$','$U(x,t_{849}=10.2)$');
set(L,'fontsize',30,'interpreter','latex')

set(gca,'fontsize',30, 'ylim',[-.5,1])








