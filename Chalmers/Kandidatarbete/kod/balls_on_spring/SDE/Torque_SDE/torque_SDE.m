%% model of two balls on a spring
clc;clf;clear;

dim=2;               %dimention
N=20;m=ones(1,N); %m(1)=100000;m(end)=10000;


%m=[1000 1 1000];N=length(m); 

%M=sum(m);            %mass of the particles

%k=[10; 100];  %spring constant
k=ones(N-1,1)*10;  %spring constant
n=0;                 %dampening
%R0=[.2; .5];              %equilibrium distance
R0=ones(N-1,1)*1;              %equilibrium distance


s     = ones(N-2,1)*20;
theta0= ones(N-2,1)*0;



%%%%%%%%%%%%%%%%%%%%% Initial condition %%%%%%%%%%%%%%%%%%%%%
%q0=[ -1,0,0,  0,.5,0, 1,0,0;     0,0,0,  0, 0,2, 0,0,0 ];
%q0=[-1.3 1, -1 0, 1 0, 1 -1.3; zeros(1,dim*N) ];
%q0=[-2 0, -1 0, 1 0, 2 0; zeros(1,dim*N) ];
r0=mean(R0);
q0=[reshape([linspace(-r0*N/2, r0*N/2, N); zeros(1,N)], 1,[]);
    zeros(1,dim*N)];



%No movement of CM
%q0=[ -.5 , 0 , -.3, .3, .1, -.3, 0, .7, 1,
%       -1 , 0 ,  0 , .7 , 1 , 0, 0, 0, 1, ];

%q0=[linspace(-1,1,N*dim);0*linspace(1,-1,N*dim)];
%q0(1,1)=0;        q0(1,2)=0;
%q0(2,1)=0;        q0(2,2)=0;
%q0(1,end-1)=-1.5; q0(1,end)=-3;
%q0(2,end-1)= 0  ; q0(2,end)=0;

%q0=randn([2,N*dim]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




dt=5e-5;
damp=mean(k)./mean(m)*dt; %this is interesting: find out more!

%T=2*pi/sqrt(k/m) %test


nVariables=N*dim;
nTrials=5/dt;




odefun=@(t, q) dq_torque(q, m, k, R0, s, theta0, damp, dim );
%odefun(0, q0)

%

brownian_amplitude=10;
G = @(t,Q) diag(repmat([0, 1], 1, nVariables))*brownian_amplitude;

SDE = sde(odefun, G, 'StartState', reshape(q0, [],1))


[Paths, Times, Z] = simByEuler(SDE, nTrials, 'DeltaTime', dt);

%
clf
disp('playing')
play_movie_v2(Times,Paths,m,dim, 20, 1)
disp('done')

%[T, Q] = ode45(odefun, [0,100], reshape(q0,1, []));
%play_movie(T,Q,m,dim, 10, 0.2)




































