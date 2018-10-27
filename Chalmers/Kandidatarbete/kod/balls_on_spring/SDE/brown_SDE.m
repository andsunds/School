%% Browninan motion SDE (BM)
clc;clear all;%clf

%bm is the constructor for the BM (brownian motion) class
obj = bm(0, 0.3) 

%dt=5e-2;
N=2;


[x, T]=obj.simByEuler(100, 'nTrials',N);%simulates 100 steps N times

plot(x(:,1),x(:,2))


%% Self constucted model
clc;clear all;clf


nVariables=1;
nTrials   =100;


X0=zeros(nVariables,1);

F = @(t,X) zeros(nVariables,1);

G = @(t,X) eye(nVariables);%diag(randn(nVariables));

SDE = sde(F, G, 'StartState', X0)


[Paths, Times, Z] = simByEuler(SDE, nTrials);


plot(Times, Paths)
%plot(Paths(:,1), Paths(:,2))



%% Ball on spring
clc;clear all;clf

nTrials =500000;

R0=0.5;
k=197;
m=17;
dt=1e-3;
damp=(k/m)*dt; %this is interesting: find out more!

%T=2*pi/sqrt(k/m) %test

Q0=[1; 0];

nVariables=length(Q0);



F = @(t,Q) ([0, 1; -k/m, -damp]*Q + [0; k*R0/m])*1;

brownian_amplitude=0;
G = @(t,Q) diag([1, 0])*brownian_amplitude;%diag(randn(nVariables));

SDE = sde(F, G, 'StartState', Q0)


% Euler forward test
% Q=zeros(nTrials,2);
% Q(1,:)=Q0';
% 
% for i=2:nTrials
%     Q(i,:)=F(i*dt,Q(i-1,:)')*dt+Q(i-1,:)';
% end
% 
% t=(0:nTrials-1)*dt;
% plot(Q(1:200,1))
% plot(t,Q(:,1))
% 
% %%

[Paths, Times, Z] = simByEuler(SDE, nTrials, 'DeltaTime', dt);


plot(Times, Paths(:,1))
%plot(Paths(:,1), Paths(:,2))










