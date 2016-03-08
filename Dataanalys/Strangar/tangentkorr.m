%% Vektoriserad tangentvektors korrelation
%<t(s) * t(s+l)> ~ exp(-l/L_P)
clc;clf;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=4;
load(['data/', filnamn{fil}])

N=size(px, 1);

n=100;%antalet punkter att kolla korr. i
K=zeros(n,1);%init.
l=linspace(0,1,n);%Vilka punkter vi ska kolla efter tangentvektor

%Beräkna tangentvektorer i alla punkter och tider. 
tangent=tangent_normal(px, py, l);

addpath('../');%Lägger till så att create_indecis kan användas
INDEX=create_indecis(n);%Tar fram index som sorterar längs med diagonalerna

tic
for i=1:N;
%                        /-- tar fram tangentvektorerna vid just tid i och
%   /-- bara övre        |    matrismultiplicerar dem för att få en matris
%   V    diagonalerna    V      med alla produkter
tmp=triu(tangent(:,:,i).'*tangent(:,:,i));

%  v-- lägger bidrag från varje bild  %v-- medelvärde över alla bilder
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






%% Tangentvektorskorr i tid 
clf;clc;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=2;
load(['data/', filnamn{fil}])

N=size(px, 1); % Total tid

n=300;%antalet punkter att kolla korr. i
K=zeros(N,n);%init.
l=linspace(0.1,.9,n);%Vilka punkter vi ska kolla efter tangentvektor

%Beräkna tangentvektorer i alla punkter och tider. 
tangent=tangent_normal(px, py, l);

for i=1:N
    T1 = tangent(:,:,i); % Tangentvektorer tid: i
    for dt=0:(N-i)
        T2 = tangent(:,:,i+dt); % Tangentvektorer tid: i+dt
        K(dt+1,:) = K(dt+1)+sum(T1.*T2, 1)./(N-dt);
        %                   ^-- summan i skalärprodukten
    end
end

Ks = sum(K,2)'/n; % Summera �ver alla punkter p� str�ngen
dt = 0:(N-1);

%figure(1)
plot(dt,(Ks))
xlabel('$d\tau$ [tid]','Interpreter','Latex');
ylabel('$<$t($\tau$)$\cdot$t($\tau+d\tau$)$>$','Interpreter','Latex')
set(gca,'Fontsize',24);%'xscale','log','yscale','log');

%%
figure(2)
surf(K(:,:)) % Summerar ej �ver punkterna p� str�ngen
xlabel('B�gl�ngd s [l�ngd]');
ylabel('dt [tid]');
zlabel('<t(s)*t(s+dt)>')
% figure(1)






