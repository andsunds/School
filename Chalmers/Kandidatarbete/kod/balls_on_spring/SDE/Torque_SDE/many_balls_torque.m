%% model of two balls on a spring
clc;clf;clear;

dim=2;               %dimension
N=4;m=ones(1,N)*.1; %m(1)=100000;m(end)=10000;

%m=[1000 1 1000];N=length(m); 

M=sum(m);            %mass of the particles

%k=[10; 100]; %spring constant
k=ones(N-1,1)*10; %spring constant
%R0=[.2; .5];%equilibrium distance
R0=ones(N-1,1)*sqrt(1); %equilibrium distance



s     = ones(N-2,1)*200;
theta0= ones(N-2,1)*0;

damp=1;%dampening

%dim==2
%q0=[ -1,0,  0,1, 1,0;  0,0, 0,0, 0,0 ];
q0=[-1.3 1, -1 0, 1 0, 1 -1.3; zeros(1,dim*N) ];

%dim==3:
%q0=[ -1,0,0,  1,1,0, 1,0,0;     0,0,0,  0, 0,0, 0,0,0 ];

%No movement of CM
%q0=[ -.5 , 0 , -.3, .3, .1, -.3, 0, .7, 1,
%       -1 , 0 ,  0 , .7 , 1 , 0, 0, 0, 1, ];

%q0=[logspace(-1,1,N*dim);linspace(1,-1,N*dim)];
%q0(1,1)=0;        q0(1,2)=0;q0(2,1)=0;        q0(2,2)=0;
%q0(1,end-1)=-1.5; q0(1,end)=-3;q0(2,end-1)= 0  ; q0(2,end)=0;

%rng(13);
%q0=randn([2,N*dim]).*[ones(1,N*dim); zeros(1,N*dim)];





odefun=@(t, q) dq_torque(q, m, k, R0, s, theta0, damp, dim );
%reshape(q0, [],8) %test
%odefun(0, q0)'%test



[T, Q] = ode45(odefun, [0,10], reshape(q0,1, []));

disp('Start movie?')
%pause
play_movie_v2(T,Q,m,dim, 8, 0.2)
%%
clc;clf
%E=total_energy(Q, m, k, dim);plot(T,E);std(E)/mean(E)

A=get_angles(Q,dim);plot(T,A/pi);grid on

































