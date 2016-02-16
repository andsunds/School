%%
clc;clf;clear all

data_energydepleted =load('energydepletedcells.csv');
%size(data_energydepleted)
data_logphase=load('logphasecells.csv');


%%
C = separera(data_logphase);

i=2; 
X=C{i}; %Examine the trajectory of the i:th particle

mean(X(:,2:3)) %Mean value of x and y coordinates
std(X(:,2:3)) %Standard deviation

figure(1)
plot(X(:,1),X(:,2:3)); %Displacement in x and y direction versus time

%%

C = separera(data_logphase);

i=5; 
X=C{i}; %Examine the trajectory of the i:th particle

figure(2)
d=diff(X(:,2:3),1,1); %Actual steps in x and y direction
hist(d) %Plot in a histogram
title('Stegfördelning i x- och y-led')

figure(3)
D=sqrt(sum(d.^2,2)); %Length of steps
hist(D)
title('Steglängdsfördelning')

%Hur medelsteget ändras med tiden för 1 partikel
i_max=length(D);
L_mean_t=zeros(1,i_max);
for i=1:i_max
    L_mean_t(i)=sum(D(1:i))/(i+1);
end
plot(X(2:end,1),L_mean_t)
%%

figure(4)
hist(data_logphase(:,4))
title('Intensitetfördelning logphase')

%%
%plot total displacement from original position for some particles
clf

figure(5)
hold on
for i=1:6; %index_largeparticle
    X=C{i};
    Intensity=X(1,4);
    %disp(i)
    %disp(Intensity)
    x_o=X(:,2)-X(1,2); %Displacement in x direction relative to original position
    y_o=X(:,3)-X(1,3); %Displacement in y direction relative to original position
    r=[sqrt(x_o.^2+y_o.^2)]; %Total displacement from original position
    r_mean=mean([sqrt(x_o.^2+y_o.^2)]);
    plot(X(:,1),r)
    %plot(i,r_mean,'*')
end
%t=1:1:10;
%plot(t,0.2*sqrt(t))
hold off
%legend('1','2','3','4','5','6','Location','Best')


%% Filmuppspelning
clf
T=X(:,1);
n=length(T);

Q=[reshape(X(:,2:3), [],1,2), zeros(n,1,2)];
m=1;
dim=2;
framesize=max(max(Q))*1.1;
playbackspeed=.05;

addpath('../balls_on_spring/universal_tools/')
play_movie_v3(T,Q,m,dim,framesize, playbackspeed,0 )


