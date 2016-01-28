%%
clc;clf;clear all

data =load('energydepletedcells.csv');
size(data);

%%
C = separera(data);

i=6;
X=C{i};

mean(X(:,2:3))
std(X(:,2:3))

plot(X(:,1),X(:,2:3));%plottar 



d=diff(X(:,2:3),1,1);%faktiska steg
hist(d)
D=sqrt(sum(d.^2,2));%steglängd
%hist(D)

%% Filmuppspelning
T=X(:,1);
n=length(T);

Q=[reshape(X(:,2:3), [],1,2), zeros(n,1,2)];
m=1;
dim=2;
framesize=max(max(Q))*1.1;
playbackspeed=.05;

play_movie_v3(T,Q,m,dim,framesize, playbackspeed,0 )


