%% Original equation of motion
clc;clf;clear


n2=.25;%.2462097;
eps=.01;

r=@(t) 1+eps*sin(t);
dr=@(t) +eps*cos(t);

%dZ=@(t,Z) [Z(2); -(+2*dr(t).*Z(2)./r(t) + n2*sin(Z(1))./r(t)) ];
dZ=@(t,Z) [Z(2); -(+2*dr(t).*Z(2)./r(t) + n2*(Z(1))./r(t)) ];
%dZ=@(t,Z) [Z(2); -(+2*dr(t).*Z(2)./r(t) + n2*(Z(1)-(Z(1).^3)/6)./r(t)) ];
    
Z0=[.1;0];
TSPAN=[0,600];

[T,X]=ode45(dZ,TSPAN,Z0);

R=r(T);
dR=dr(T);
Y=X(:,1);
dY=X(:,2);
%plot(R.*sin(Y),-R.*cos(Y))

%
plot(T,Y)
hold on
plot(T,Z0(1)*exp(3*eps*T/8))
%plot(T,cos(T),'--')
%{
This does not look very promisin. The approximate version seems to result
in much more well behaved solutions.
%}

%% Energy
clf;clc
%E=.5*(dR.^2+(R.*dY).^2)+n2*(1-R.*cos(Y));
E=.5*(0*dR.^2+(R.*dY).^2)+n2*R.*(1-cos(Y));

plot(T,log(E))
hold on
plot(T,log(E(1))+(3/4)*eps*T)
set(gca,'ylim',[-7,-2])

%% Saving simulations
data=[T,Y,E];
save('swing_sim_optimal_lin.tsv','data','-ascii', '-double', '-tabs')



%% No simulation, just a plot of how we can choose r(t)
clc;clf;clear;

phi0=10*pi/180;

phi=@(t) phi0*cos(t);
dphi=@(t) -phi0*sin(t);

eps=.1;
r=@(x,v) 1+eps*((x-v).^2./(x.^2+v.^2));
%r=@(x,v) 1+eps*cos(2*atan2(x,v)+pi/4+0*pi/180);



x=@(t) r(phi(t),dphi(t)).*sin(phi(t));
y=@(t) -r(phi(t),dphi(t)).*cos(phi(t));

%{
t=linspace(0,2*pi-0.3,1000);
plot(x(t),y(t)), hold on
plot(sin(phi(t)),-cos(phi(t)),'k--')
axis([-.25,.25,-1.1,-.9])
axis equal
%}
% {
t=linspace(0,4*pi,1000);
T=t/(2*pi);
R=r(phi(t),dphi(t));
plot(T,(R-mean(R))/max(R))
hold on
plot(T,phi(t))
plot(T(1:end-1), diff(R)*50,'--')

grid on
%}


%% Forcing depending on the angle
%% No good
%Stability test for different a's and b's
clc;clf;clear

phi0=10*pi/180;%inital angle
eps=.1;

A=(0.1:0.025:1); %angle factor
B=(-180:2.5:180); %phase lag

ampl_ratio=zeros(length(B),length(A));

%loop over different combinations of a and b
for i=1:length(A)
    a=A(i)*pi/phi0;
    for j=1:length(B)
        b=B(j)*pi/180;
        %cos driving
        r=@(t) 1+eps*cos(a*t+b);%radius of swing
        dr=@(t) -a*eps*sin(a*t+b);%t derivative of the radius

        %sin driving
        %r=@(t) 1+eps*sin(a*t).^2;
        %dr=@(t) +a*eps*sin(2*a*t);

        dZ=@(t,Z) [Z(2); -2*dr(Z(1))./r(Z(1))*(Z(2)).^2 - 1/r(Z(1)).*Z(1)];
    
        Z0=[phi0;0];
        TSPAN=[0,500];
        [T,X]=ode45(dZ,TSPAN,Z0);
        
        L_T=length(T);
        
        ampl_ratio(j,i)=max(abs(X(floor(0.9*L_T):end,1)))/phi0;
    end
end

%%
%A=(0.1:0.025:1);
%B=(-30:1.25:30);
Ia=find(A<=2/pi);
Ib=find((B<=180)&(B>=-180));
surf(A(Ia),B(Ib),log(ampl_ratio(Ib,Ia)))
ylabel('b')
xlabel('a')

%%

Y=X(:,1);%/Z0(1);
R=r(X(:,1));

%figure(1);clf
plot(T,Y*180/pi), hold on

%plot(X(:,1),X(:,2));axis equal

%% Plot of the CM's movement
clf

plot(R.*sin(Y),-R.*cos(Y))
hold on

plot([-1,1], [0,0],'k')

axis equal







