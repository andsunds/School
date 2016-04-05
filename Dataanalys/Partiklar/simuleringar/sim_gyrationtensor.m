%% Simuleringar
clc;clf; clear all

N_steps=1000;
N_trials=100000;
tic
for skev=[1]
kvot=zeros(N_trials,1);
egen=zeros(N_trials,2);
for i=1:N_trials
    X=(rand(N_steps,1)*2-1);
    Y=(rand(N_steps,1)*2-1)*skev;
    gyr = [sum(X.^2), sum(X.*Y); 
           sum(Y.*X), sum(Y.^2)]/N_steps;
       
    %gyr = [sum(X(:,i).^2), sum(X(:,i).*Y(:,i)); sum(Y(:,i).*X(:,i)), sum(Y(:,i).^2)]/N_steps;
    [~,D] = eig(gyr);

    kvot(i)=(D(4)/D(1))^((randn>0)*2 -1 );
    egen(i,:)=[D(1), D(4)];
end
%sparar sim
%{
if skev==1
    save('kvot.mat', 'kvot', 'egen', '-mat')
else
    save(sprintf('kvot_skev_%d.mat', skev), 'kvot', 'egen', '-mat')
end
%}
end
toc
%%
clf,clc
[N_X, X]=hist(kvot,100);
hold on
plot(X+0*.5*mean(diff(X)), (N_X)/sum(N_X), '.')

G=@(n) (exp(-1)/sqrt(2*pi)) * (2*(n-1)/(n-2))^(n-1)*sqrt(n-1);

f=@(w, n) G(n) .*(sqrt(w)./(1+w)).^n ./w;
            %.* w.^(n/2-1).*(1+w).^(-n);
x=linspace(0,1.5);

plot(x,f(x,1000)/200)



%clear all
%%
clc;clear all;clf


b=dir('kvot*.mat');
leg_str=cell(length(b),1);
n_bins=1000;
for i=1:length(b)
    leg_str{i}=b(i).name;
    load(leg_str{i})
    [N_X, X]=hist(kvot,n_bins);
    hold on
    plot(X+.5*mean(diff(X)), cumsum(N_X)/sum(N_X), '.')
    %plot(X, N_X/sum(N_X), '*');hold on
end
l=legend(leg_str);set(l, 'fontsize', 20)
set(gca, 'fontsize', 20, 'yscale', 'log', 'xscale', 'lin')


