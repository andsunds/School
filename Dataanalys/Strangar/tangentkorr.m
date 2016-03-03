%% Tangentvektorskorr (Finns förbättringspotential)
%<t(s) * t(s+l)> ~ exp(-l/L_P)
clf;clc;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=3;
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
%% Vektorisering
clc;clf;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=3;
load(['data/', filnamn{fil}])

N=size(px, 1);

n=100;%antalet punkter att kolla korr. i
K=zeros(n,1);%init.
l=linspace(0,1,n);%Vilka punkter vi ska kolla efter tangentvektor

%Beräkna tangentvektorer i alla punkter och tider. 
tangent=tangent_normal(px, py, l);

addpath('../Partiklar/');%Lägger till så att create_indecis kan användas
INDEX=create_indecis(n);%Tar fram index som sorterar längs med diagonalerna
tic
for i=1:N;
%                        /-- tar fram tangentvektorerna vid just tid i och
%   /-- bara övre        |    matrismultiplicerar dem för att få en matris
%   V    diagonalerna    V      med alla produkter
tmp=triu(tangent(:,:,i).'*tangent(:,:,i));

%  v-- lägger bidrag från varje bild  %v-- medelvärde över alla bilder
K=K+sum( tmp(INDEX), 2)./fliplr(1:n).'/N;
%                       ^normerar m.a.p. antalet datapunkter
end
toc



plot(linspace(0,1,n), K)

















