%% Kovariansundersökning
clc;clf;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=1;
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
Q0=bsxfun(@minus, Q0, mean([polyval(PX_mean, S); polyval(PY_mean, S)], 2) );%tyngdpunkt i origo

for i=1:N %loopa över alla bilder (tid)
    Q1=[polyval(px(i,:), P); polyval(py(i,:), P)]; 
    Q1=bsxfun(@minus, Q1, mean([polyval(px(i,:), S); polyval(py(i,:), S)], 2) );%tyngdpunkt i origo
    
    A(:,i)  = diag(N0.'*(Q1-Q0)); %tidsutv. i dim. 2
    M=M+A(:,i)*A(:,i).';
end
toc
%A_mean=sum(A, 2)/n;%medelvärde av summan
M=M/N;


tic
[V, D] = eig(M);
toc
size(V)

B=V.'*A;
size(B)


plot(B.')
