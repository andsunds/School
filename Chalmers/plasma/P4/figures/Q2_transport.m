%%
clc;clear
load('constants.mat');

sigma=1e-19;
nN=64*pi*epsilon_0^2*sigma/e^2

N=1e26;
D=1/(N*sigma)*sqrt(2*e/m_e)

%%
clc;clf;clear

a1=0.05; a2=0.5; a3=5;

S0=10; r0=0.1; R=1;
T0=1e-1;


chi1=@(T) a1./sqrt(T);
r1=r0/2;
chi2=@(r,T) a2./sqrt(T)*r^(-3/2),%(r1^(-3/2)+(r.^(-3/2)-r1^(-3/2)).*heaviside(r-r1));
chi3=@(T) a3*(T.^(3/2));

[RR1,TT1]=ode45(@(r,T) -S0*heaviside(r-r0)./chi1(T),[R,0],T0);
[RR2,TT2]=ode45(@(r,T) -S0*heaviside(r-r0)./chi2(r,T),[R,0],T0);
[RR3,TT3]=ode45(@(r,T) -S0*heaviside(r-r0)./chi3(T),[R,0],T0);

plot(RR1,TT1,'--',RR2,TT2,'--',RR3,TT3,'--','linewidth',3), hold on

r=linspace(0,1);
t1=(sqrt(T0)+S0/(2*a1)*((R-r0)-(r-r0).*heaviside(r-r0))).^2;
t2=(sqrt(T0)+S0/(5*a2)*((R^(5/2)-r0^(5/2))-(r.^(5/2)-r0.^(5/2)).*heaviside(r-r0))).^2;
t3=(T0^(5/2)+5*S0/(2*a3)*((R-r0)-(r-r0).*heaviside(r-r0))).^(2/5);
plot(r,t1,r,t2,r,t3,'linewidth',1)

t1(1),t2(1),t3(1),

FS=20;
set(gca,'fontsize',FS,'yscale','log','xscale','lin','ylim',[1e-1,1e4])

plot([T0,T0],[2.2e-1,1e5],'k:')
text(0.8,6.e2,'Classical','fontsize',FS,'interpreter','latex')
text(0.65,1.3e1,'Neoclassical','fontsize',FS,'interpreter','latex')
text(0.5,2.2e0,'Gyro-Bohm','fontsize',FS,'interpreter','latex')
text(0.09,1.5e-1,'$r_0$','fontsize',FS,'interpreter','latex')
xlabel('$r$','fontsize',FS,'interpreter','latex')
ylabel('$T(r)$','fontsize',FS,'interpreter','latex')