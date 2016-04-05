%% Find noise
% x_mät = x + brus
clc;clf;clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=1;
C=separera(load(filnamn{fil}));

n=1000;
dt=10e-3;
f=(0:(n/2-1))'/dt/(n-1);

Spektr=zeros(n,2);

for i=find(cellfun('length',C)>=n).'; 
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    %TN=C{i}(:,2:3);%laddar in data för partikeln
    
    Spektr=Spektr+abs(fft(TN, [],1)).^2;
end

S=sum(Spektr(1:(n/2),:), 2);

plot(f,S), hold on
xlabel('$f/[\mathrm{Hz}]$','interpreter', 'latex')
ylabel('Spektrum','interpreter', 'latex')
set(gca,'FontSize',15,'XScale','log','YScale','log')

show=101;
c=[log(f(2:show)), ones(show-1,1)]\log(S(2:show));

x=logspace(-1,2);
c(1)
plot(x, exp(c(2))*x.^c(1))

%%
clc;clf;clear all

n=1024;
dt=10e-3;
f=(0:(n/2-1))'/dt/(n-1);

n_sim=1000;

X=3e-8*randn(n,2*n_sim);%laddar in data för partikeln

Spektr=abs(fft(X, [],1)).^2;


S=sum(Spektr(1:(n/2),:), 2)/n_sim;

plot(f,S), hold on
xlabel('$f/[\mathrm{Hz}]$','interpreter', 'latex')
ylabel('Spektrum','interpreter', 'latex')
set(gca,'FontSize',15,'XScale','lin','YScale','log')

mean(S)















