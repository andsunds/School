%% Hitta hur stort bruset är i datan, samt plottar 
% x_mät = x + brus
clc;clf;clear all
addpath('../');

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=2;
C=separera(load( filnamn{fil} ));

n=992;
dt=10e-3;
f=(0:(n/2-1))'/dt/(n-1);

Spektr=zeros(n,2);
VAR=0;

index=find(cellfun('length',C)>=n).'; 
tic
for i=index;
    TN=koordinatbyte(C{i}(1:n,2:3));%laddar in data för partikeln
    %TN=C{i}(:,2:3);%laddar in data för partikeln
    
    VAR=VAR+sum(var(diff(TN,1,1),0,1),2);
    
    Spektr=Spektr+abs(fft(TN, [],1)).^2;
end
toc
S=sum(Spektr(1:(n/2),:), 2)/length(index)/(n);

steglangd=sqrt(VAR/length(index));

%S=S-(11e-9)^2;

%save('spektrum_logphase.mat', 'S', 'f', '-mat')



fprintf('Standardavvikelsen bland stegen var i medel: %3.3f nm\n\n',steglangd*1e9);

plot(f,S), hold on
xlabel('$f/[\mathrm{Hz}]$','interpreter', 'latex')
ylabel('Spektrum /[m$^2$/Hz]','interpreter', 'latex')
set(gca,'FontSize',15,'XScale','log','YScale','log')
grid on

show=101;
c=[log(f(2:show)), ones(show-1,1)]\log(S(2:show));

x=logspace(-1,1.6);
fprintf('Anpassad lutning till de första %2.0f %% av spektrat var: %2.2f\n\n',...
         100*show/n, c(1))
plot(x, exp(c(2))*x.^c(1))




%% Simulera Wienerprocess med olika muckat brus
clc;clf;clear all

n=2^10;
dt=10e-3;
f=(0:(n/2-1))'/dt/(n-1);

n_sim=1000;

sigma_brus=.5e-9;

X=2.256e-9*cumsum(randn(n,n_sim)) + sigma_brus*randn(n,n_sim);%laddar in data för partikeln

Spektr=abs(fft(X, [],1)).^2;


S=sum(Spektr(1:(n/2),:), 2)/n_sim/n;

%S=S-(pi*sigma_brus)^2;

plot(f,S), hold on
xlabel('$f/[\mathrm{Hz}]$','interpreter', 'latex')
ylabel('Spektrum','interpreter', 'latex')
set(gca,'FontSize',15,'XScale','log','YScale','log')
%axis([1e-1, 1e2, 1e-18, 1e-12])
grid on


%% simulera O-U
clc;clf;clear all


N_steps=1024;
N_trials=10000;

dt=10e-3;
f=(0:(N_steps/2-1))'/dt/(N_steps-1);

Spektr=zeros(N_steps,2);

tau=[2.6, 6];

for l=1:length(tau)
tic
k=1e-2./tau(l);

for i=1:N_trials
    W=2.256e-9*randn(N_steps-1,2);
    X=zeros(N_steps,1);
    Y=zeros(N_steps,1);
    for j=2:N_steps
        X(j)=(1-k)*X(j-1)+W(j-1,1);
        Y(j)=(1-k)*Y(j-1)+W(j-1,2);
    end
    Spektr=Spektr+abs(fft([X,Y], [],1)).^2;
    
end

S=sum(Spektr(1:(N_steps/2),:), 2)/N_trials/N_steps;
toc

% % save(sprintf('spektrum_%d_spring_%1.1f.mat', N_steps, tau(l)) , 'S', 'f');



hold on
plot(f,S)
set(gca,'FontSize',15,'XScale','log','YScale','log')
pause(.1)
end

%% Plottar ALLT
clc;clf;clear all;
hold on

b=dir('spektrum*.mat');
leg_str=cell(length(b),1);

for i=1:length(b)
    leg_str{i}=b(i).name;
    load(leg_str{i})
    plot(f,S)
    
    index=find(f>=2e-1 & f<=1e1);
    c=[log(f(index)), ones(length(index),1)]\log(S(index));
    fprintf('%s, lutning: \t%1.2f \n', leg_str{i}, c(1))
end

axis([2e-3,5e1, 1e-18, 1e-11])

l=legend(leg_str);set(l, 'fontsize', 20)
xlabel('$f/[\mathrm{Hz}]$','interpreter', 'latex')
ylabel('Spektrum /[m$^2$/Hz]','interpreter', 'latex')
set(gca,'FontSize',15,'XScale','log','YScale','log')
grid on
















