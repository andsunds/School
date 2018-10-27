%%
clc;clear,clf

NUn=[0,2];
NUT=linspace(0,4,128);

J=length(NUn);
K=length(NUT);

Tmin=zeros(K,J);
Pmin=zeros(K,J);

for j=1:J % for every \nu_n
    nun=NUn(j);
    parfor k=1:K % for every \nu_T
        %Find at what average temp \ev{T}_V we get minimum value of P\tau_E
        [Tmin(k, j),Pmin(k,j)]=fminsearch(@(TV)RHS_ignition(TV, NUT(k), nun), 5);
    end
end
%renormalizing to get actual P\tau_E
KL=0.15; Ka=1.37;
Pmin=Pmin*KL/Ka;
%% Plotting
clc;clf
plot(NUT,Tmin,'linewidth',1.5)
set(gca,'ylim',[0,15],'fontsize',18)
xlabel('$\nu_T$','interpreter','latex')
ylabel('$T_\mathrm{min}$','interpreter','latex')

l=legend('$\nu_n=0$','$\nu_n=2$');
set(l,'interpreter','latex')

%%
clc;clf;clear

KL=0.15; Ka=1.37;
NUT=linspace(2,4,4);
nun=0;
TV=linspace(0,20);
J=length(TV);
K=length(NUT);

P=zeros(J,K);


parfor j=1:J
for k=1:K
    P(j,k)=RHS_ignition(TV(j), NUT(k), nun)*KL/Ka;
end
end

plot(TV,P)
set(gca,'yscale','log')



%%
clc;clear;clf

x=linspace(0,1);
nu=40;
y=(1+nu)*(1-x.^2).^(nu);

plot(x,y.*x)





