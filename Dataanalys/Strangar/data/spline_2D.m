%% Detta verkar ge en spline kurva.
clc;clf;clear all
load('nonconfined_167min.mat')
i=1;
t=1:N_points(i);
xy = [coordinates(i,t,1); coordinates(i,t,2)];

sp=csape(t,xy,'spline');%vet inte exat vad den här gör men det blir mjukt i alla fall

%fnplt(sp,1)

XY=ppval(linspace(0,365,10000), sp);
%size(XY)

plot(XY(1,:), XY(2,:))
