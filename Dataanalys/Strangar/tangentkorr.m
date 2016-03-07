%% Vektoriserad tangentvektors korrelation
%<t(s) * t(s+l)> ~ exp(-l/L_P)
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

%Ber√§kna tangentvektorer i alla punkter och tider. 
tangent=tangent_normal(px, py, l);

addpath('../');%L√§gger till s√• att create_indecis kan anv√§ndas
INDEX=create_indecis(n);%Tar fram index som sorterar l√§ngs med diagonalerna

tic
for i=1:N;
%                        /-- tar fram tangentvektorerna vid just tid i och
%   /-- bara √∂vre        |    matrismultiplicerar dem f√∂r att f√• en matris
%   V    diagonalerna    V      med alla produkter
tmp=triu(tangent(:,:,i).'*tangent(:,:,i));

%  v-- l√§gger bidrag fr√•n varje bild  %v-- medelv√§rde √∂ver alla bilder
K=K+sum( tmp(INDEX), 2);
%                       ^normerar m.a.p. antalet datapunkter
end
K=K./fliplr(1:n).'/N;
toc


plot(linspace(0,1,n), K)




%% Tangentvektorskorr bara for-loopar
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

%Ber√§kna tangentvektorer i alla punkter och tider. 
tangent=tangent_normal(px, py, l);

tic
for i=1:N;
T=tangent(:,:,i);%tar fram tangentvektorerna vid just tid i.

%Ber√§kna korr
for s=1:n
for dl=1:(n-s)
    %Sumerar delar till medelv√§rden av korrelationsfunktionen
    K(dl)=K(dl)+(T(:,s).'*T(:,s+dl-1))/(n-dl)/N;
end
end

end
toc
%plottar
plot(linspace(0,1,n-1), K)

%axis([0,1, 0,1])

%set(gca, 'fontsize',15, 'yscale', 'log')






%% Tangentvektorskorr i tid 
clf;clc;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=4;
load(['data/', filnamn{fil}])

N=size(px, 1); % Total tid

n=100;%antalet punkter att kolla korr. i
K=zeros(N,n);%init.
l=linspace(0,1,n);%Vilka punkter vi ska kolla efter tangentvektor

%Ber√§kna tangentvektorer i alla punkter och tider. 
tangent=tangent_normal(px, py, l);

for i=1:N
    T1 = tangent(:,:,i); % Tangentvektorer tid i
    for dt=0:(N-i)
        T2 = tangent(:,:,i+dt); % Tangentvektorer tid i+dt
        K(dt+1,:) = K(dt+1)+sum(T1.*T2)./(N-dt);
    end
end

Ks = sum(K,2)'/n; % Summera ˆver alla punkter pÂ str‰ngen
dt = linspace(0,N-1,N);
figure(1)
plot(dt,Ks)
xlabel('$d\tau$ [tid]','Interpreter','Latex');
ylabel('$<$t($\tau$)$\cdot$t($\tau+d\tau$)$>$','Interpreter','Latex')
set(gca,'Fontsize',24);%'xscale','log','yscale','log');


figure(2)
surf(K(:,:)) % Summerar ej ˆver punkterna pÂ str‰ngen
xlabel('BÂgl‰ngd s [l‰ngd]');
ylabel('dt [tid]');
zlabel('<t(s)*t(s+dt)>')
% figure(1)






