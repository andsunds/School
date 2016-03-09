%% Kovariansundersökning
clc;clf;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=4;
load(['data/', filnamn{fil}])

N=size(px, 1);
S=linspace(0,1,1000);

n=100;%antalet punkter att kolla korr. i

A=zeros(n,N);%init
M=zeros(n);

tic
P=linspace(0,1,n);%loopa över hela strängen

%Beräkna normalvektorer i alla punkter och tider. 
[~, N0]=tangent_normal(PX_mean, PY_mean, P);

Q0=[polyval(PX_mean, P); polyval(PY_mean, P)];%Den undersökt punkten på medelkurvan
Q0=bsxfun(@minus, Q0,...
    mean([polyval(PX_mean, S); polyval(PY_mean, S)], 2) );%tyngdpunkt i origo

for i=1:N %loopa över alla bilder (tid)
    Q1=[polyval(px(i,:), P); polyval(py(i,:), P)]; 
    Q1=bsxfun(@minus, Q1,... 
        mean([polyval(px(i,:), S); polyval(py(i,:), S)], 2) );%tyngdpunkt i origo
    
    A(:,i)  = diag(N0.'*(Q1-Q0)); %tidsutv. i dim. 2
    M=M+A(:,i)*A(:,i).';
end
M=M/N;%medelvärde
toc

tic
[V, ~] = eig(M);%tar fram egenvektorer
toc


B=V.'*A;
%B = [B1(1) B1(2), B1(3) ...
%     B2(1) B1(2), B1(3) ...
%     B3(1) B3(2), B3(3) ...
%     ...                   ]
% Där Bi(t) är den i:te moden i bild/tid t.


figure(1)
plot(A.') %Massa bröte bara
figure(2)
plot(B.') %mindre bröte, men ändå otydligt

%% Autokorrelation
K=zeros(N,n);%init.
%Lägger till så att create_indecis kan användas
addpath('../');
%Tar fram index som sorterar längs med diagonalerna
INDEX=create_indecis(N);

tic
for i=1:n;
%   /-- bara övre        
%   V    diagonalerna
tmp=triu(B(i,:).'*B(i,:));
%                ^-- korrelationsfunktion på diagonalerna
%                         %v-- medelvärde över diagonalerna
K(:,i)=sum( tmp(INDEX), 2)./fliplr(1:N).';
end
toc


figure(3)
plot(0:(N-1), K) %Häftig korrelatinosfunktioner

%% Korskorrelation
clc;clf
K=zeros(N,1);%init.
%Lägger till så att create_indecis kan användas
addpath('../');
%Tar fram index som sorterar längs med diagonalerna
INDEX=create_indecis(N);

i=90;
j=90;

tic
%   /-- bara övre        
%   V    diagonalerna
tmp=triu(B(i,:).'*B(j,:));
%                ^-- korrelationsfunktion på diagonalerna
%                         %v-- medelvärde över diagonalerna
K_kors=sum( tmp(INDEX), 2)./fliplr(1:N).';

tmp=triu(B(i,:).'*B(i,:));
K_auto1=sum( tmp(INDEX), 2)./fliplr(1:N).';

tmp=triu(B(j,:).'*B(j,:));
K_auto2=sum( tmp(INDEX), 2)./fliplr(1:N).';
toc


%figure(3)
plot(0:(N-1), K_kors) 
hold on
plot(0:(N-1), K_auto1) 
plot(0:(N-1), K_auto2) 


