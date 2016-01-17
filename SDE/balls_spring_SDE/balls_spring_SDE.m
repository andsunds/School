%% model of two balls on a spring
clc;clf;clear;

dim=2;               %dimention
N=2;m=ones(1,N); %m(1)=100000;m(end)=10000;


%m=[1000 1 1000];N=length(m); 

%M=sum(m);            %mass of the particles

%k=[10; 100];  %spring constant
k=ones(N-1,1)*100;  %spring constant
n=0;                 %dampening
%R0=[.2; .5];              %equilibrium distance
R0=ones(N-1,1)*.1;              %equilibrium distance



%%%%%%%%%%%%%%%%%%%%% Initial condition %%%%%%%%%%%%%%%%%%%%%
%q0=[ -1,0,0,  0,.5,0, 1,0,0;     0,0,0,  0, 0,2, 0,0,0 ];

%No movement of CM
%q0=[ -.5 , 0 , -.3, .3, .1, -.3, 0, .7, 1,
%       -1 , 0 ,  0 , .7 , 1 , 0, 0, 0, 1, ];

q0=[linspace(-1,1,N*dim);0*linspace(1,-1,N*dim)];
%q0(1,1)=0;        q0(1,2)=0;
%q0(2,1)=0;        q0(2,2)=0;
%q0(1,end-1)=-1.5; q0(1,end)=-3;
%q0(2,end-1)= 0  ; q0(2,end)=0;

%q0=randn([2,N*dim]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




dt=1e-3;
damp=(k./m)*dt %this is interesting: find out more!

%T=2*pi/sqrt(k/m) %test


nVariables=N*dim;
nTrials=10000;
dt=1e-3;
damp=(k./m)*dt; %this is interesting: find out more!


odefun=@(t, q) dq_SDE(q, m, k, damp, R0, N, dim);
odefun(0, q0)

%%

brownian_amplitude=0;
G = @(t,Q) diag(repmat([1, 0], 1, nVariables))*brownian_amplitude;

SDE = sde(odefun, G, 'StartState', reshape(q0, [],1))


[Paths, Times, Z] = simByEuler(SDE, nTrials, 'DeltaTime', dt);

play_movie(Times,Paths,m,dim, 10, 1)

%[T, Q] = ode45(odefun, [0,100], reshape(q0,1, []));
%play_movie(T,Q,m,dim, 10, 0.2)




































