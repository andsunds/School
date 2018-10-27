%% model of two balls on a spring
clc;clf;clear;

dim=3;               %dimention
N=50;m=ones(1,N); %m(1)=100000;m(end)=10000;

%m=[1000 1 1000];N=length(m); 

M=sum(m);            %mass of the particles

%k=[10; 100];  %spring constant
k=ones(N-1,1)*100;  %spring constant
n=0;                 %dampening
%R0=[.2; .5];              %equilibrium distance
R0=ones(N-1,1)*.1;              %equilibrium distance



%q0=[ -1,0,0,  0,.5,0, 1,0,0;     0,0,0,  0, 0,2, 0,0,0 ];

%No movement of CM
%q0=[ -.5 , 0 , -.3, .3, .1, -.3, 0, .7, 1,
%       -1 , 0 ,  0 , .7 , 1 , 0, 0, 0, 1, ];

%q0=[logspace(-1,1,N*dim);linspace(1,-1,N*dim)];
%q0(1,1)=0;        q0(1,2)=0;q0(2,1)=0;        q0(2,2)=0;
%q0(1,end-1)=-1.5; q0(1,end)=-3;q0(2,end-1)= 0  ; q0(2,end)=0;

q0=randn([2,N*dim]);





odefun=@(t, q) dq_3D(q, m, k, n, R0, N, dim);
%reshape(q0, [],8) %test
%odefun(0, q0)'%test



[T, Q] = ode45(odefun, [0,100], reshape(q0,1, []));

clf
play_movie_v2(T,Q,m,dim, 1, 0.2)

%E=total_energy(Q, m, k, dim);plot(T,E);std(E)/mean(E)



































