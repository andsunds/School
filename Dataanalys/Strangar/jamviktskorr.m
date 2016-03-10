%% Kollar korrelatin i jämvikten (bygger på koden i jamvikt2)
% korr(ds)=<(n0(s)*a(s))(n0(s+ds)*a(s+ds))>
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

Kn=zeros(n,1);%init

addpath('../');%Lägger till så att create_indecis kan användas
INDEX=create_indecis(n);%Tar fram index som sorterar längs med diagonalerna



punkter=linspace(0.1,.9,n);%loopa över hela strängen

%Beräkna normalvektorer i alla punkter och tider. 
[~, N0]=tangent_normal(PX_mean, PY_mean, punkter);

Q0=[polyval(PX_mean, punkter); polyval(PY_mean, punkter)];%medelkurvan
Q0=Q0-mean([polyval(PX_mean, S); polyval(PY_mean, S)], 2);%tyngdpunkt i origo

tic
for i=1:N %loopa över alla bilder
    Q1=[polyval(px(i,:), punkter); polyval(py(i,:), punkter)]; 
    Q1=Q1-mean([polyval(px(i,:), S); polyval(py(i,:), S)], 2);%tyngdpunkt i origo
    
    avst=(Q1-Q0);%avstondsvektor
    %avst=avst./repmat(sqrt(sum(avst.^2, 1)), 2,1);%normering
    
    normalavst=diag(N0.'*avst);
    %    v-- fancy korrelationsberäknig
    tmp=triu(normalavst*normalavst.');
    Kn=Kn+sum( tmp(INDEX), 2);
end
Kn=Kn./fliplr(1:n).'/N;
toc

plot(punkter,Kn)

xlabel('$\Delta l$ /[px]', 'interpreter', 'LaTeX')
ylabel('$<(\mathbf{r}-\mathbf{r}_0)\cdot\mathbf{\hat{n}}_0>_t$ /[px]', 'interpreter', 'LaTeX')

set(gca, 'fontsize', 20)