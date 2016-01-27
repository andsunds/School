clc;clf;clear;

dim=3;               %dimention
%N=3;

m=[1 4 .5];
%m=ones(1,N);
%m(1)=100000;m(end)=10000;

M=sum(m);            %mass of the particles
N=length(m);         %number of particles
k=[10; 100];  %spring constant
%k=ones(N-1,1)*100;  %spring constant
n=0;                 %dampening
R0=[.2; .5];              %equilibrium distance
%R0=ones(N-1,1)*.1;              %equilibrium distance
Br=0;                %Browninan impact strength



%No movement of CM
q0=[ -.5 , 0 , -.3, .3, .1, -.3, 0, .7, 1
     -1 , 0 ,  0 , .7 , 1 , 0, 0, 0, 1, ];
%%%%%%%%%%%%%%%%

dt=1e-3;

damp=(k./m)*dt %this is interesting: find out more!


q=reshape(q0, 2,[]);%makes a matrix from the vector

r=reshape(q(1,:), dim,[])';


f=springforce_SDE( r,  k, R0 );


a=eye(N,N-1)-[zeros(1,N-1); eye(N-1)];

F=a*f;

A=F./repmat(m, dim, 1).';

size(q(2,:))
size(repmat(m, dim, 1))

Dq=[q(2,:);
     reshape(A.', 1, [])-damp.*q(2,:)];

Dq=reshape(Dq, [], 1);













