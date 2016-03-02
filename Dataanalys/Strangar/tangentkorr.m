
%% Tangentvektorskorr (Finns förbättringspotential)
%<t(s) * t(s+l)> ~ exp(-l/L_P)
clf;clc;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=4;
load(['data/', filnamn{fil}])

N=size(px, 1);

n=100;%antalet punkter att kolla korr. i
K=zeros(1,n-1);%init.
l=linspace(0,1,n);%Vilka punkter vi ska kolla efter tangentvektor

%Beräkna tangentvektorer i alla punkter och tider. 
tangent=tangent_normal(px, py, l);

tic
for i=1:N;
T=tangent(:,:,i);%tar fram tangentvektorerna vid just tid i.

%Beräkna korr
for s=1:n
for dl=1:(n-s)
    %Sumerar delar till medelvärden av korrelationsfunktionen
    K(dl)=K(dl)+(T(:,s).'*T(:,s+dl-1))/(n-dl)/N;
end
end

end
toc
%plottar
plot(linspace(0,1,n-1), K)

%axis([0,1, 0,1])

%set(gca, 'fontsize',15, 'yscale', 'log')
