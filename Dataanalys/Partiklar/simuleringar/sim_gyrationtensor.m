%% Simuleringar, Wienerprocess
clc;clf; clearvars
hold on
N_steps=1000;
N_trials=100000;

%for brus=[9, 14];
brus=0;
tic

kvot=zeros(N_trials,1);
egen=zeros(N_trials,2);
for i=1:N_trials
    X=cumsum(randn(N_steps,1)) + brus*randn(N_steps,1);
    Y=cumsum(randn(N_steps,1)) + brus*randn(N_steps,1);
    gyr = [sum(X.^2), sum(X.*Y); 
           sum(Y.*X), sum(Y.^2)]/N_steps;
       
    %gyr = [sum(X(:,i).^2), sum(X(:,i).*Y(:,i)); sum(Y(:,i).*X(:,i)), sum(Y(:,i).^2)]/N_steps;
    [~,D] = eig(gyr);

    kvot(i)=( (D(4)-D(1))/(D(4)+D(1)) )^2;%^((randn>0)*2 -1 );
    egen(i,:)=[D(1), D(4)];
end
%sparar sim
% % save(sprintf('kvot_brus_%2.1f.mat', brus), 'kvot', 'egen', '-mat')

toc
%

bins=200;%logspace(0,2);
[N_X, X]=hist(kvot,bins);
plot([0, X+.5*[diff(X),0]], 1- [0, cumsum(N_X)]/sum(N_X), '-')
set(gca, 'fontsize', 20, 'yscale', 'lin', 'xscale', 'lin')
%end

%%
leg_str=cell(2,1);
leg_str{1}='kvot_energydepleted.mat';
leg_str{2}='kvot_logphase.mat';


for i=1:2
    load(leg_str{i})
    [N_X, X]=hist(kvot,bins);
    hold on
    plot(X+.5*[diff(X),0], 1-cumsum(N_X)/sum(N_X), '-o')
    %plot(X, N_X/sum(N_X), '*');hold on
end

axis([1,1e2, 1e-3, 1])

set(gca, 'fontsize', 20, 'yscale', 'log', 'xscale', 'log')



%% Simuleringar, Ornstein-Uhlenbeck
clc;clf; clearvars

N_steps=1000;
N_trials=1000000;

tau=[2.6];

for l=1:length(tau)
tic
k=1e-2./tau(l);
kvot=zeros(N_trials,1);
egen=zeros(N_trials,2);
for i=1:N_trials
    W=randn(N_steps-1,2);
    X=zeros(N_steps,1);
    Y=zeros(N_steps,1);
    for j=2:N_steps
        X(j)=(1-k)*X(j-1)+W(j-1,1);
        Y(j)=(1-k)*Y(j-1)+W(j-1,2);
    end
    gyr = [sum(X.^2), sum(X.*Y); 
           sum(Y.*X), sum(Y.^2)]/N_steps;
       
    %gyr = [sum(X(:,i).^2), sum(X(:,i).*Y(:,i)); sum(Y(:,i).*X(:,i)), sum(Y(:,i).^2)]/N_steps;
    [~,D] = eig(gyr);

    kvot(i)=(D(4)/D(1));
    egen(i,:)=[D(1), D(4)];
end
toc

bins=logspace(0,2);
[N_X, X]=hist(kvot,bins);
hold on
plot(X+.5*[diff(X),0], 1-cumsum(N_X)/sum(N_X), '-*')
pause(0.1)
% % save(sprintf('kvot_spring_%1.1fs.mat', tau(l)), 'kvot', 'egen', '-mat')

end
set(gca, 'fontsize', 20, 'yscale', 'log', 'xscale', 'log')
%

leg_str=cell(3,1);
leg_str{1}='kvot_energydepleted.mat';
leg_str{2}='kvot_logphase.mat';
leg_str{3}='kvot.mat';

for i=1:2
    load(leg_str{i})
    [N_X, X]=hist(kvot,bins);
    hold on
    plot(X+.5*[diff(X),0], 1-cumsum(N_X)/sum(N_X), '-o')
    %plot(X, N_X/sum(N_X), '*');hold on
end

axis([1,1e2, 1e-3, 1])

set(gca, 'fontsize', 20, 'yscale', 'log', 'xscale', 'log')

%% Plottar ALLT, egenvärdeskvot
clc;clearvars;clf


b=dir('kvot*.mat');
leg_str=cell(length(b),1);
L=length(b);
bins=logspace(0,2,50);

data=zeros(length(bins)+1, 2*L);


for i=1:L
    leg_str{i}=b(i).name;
    load(leg_str{i}, 'kvot')
    [N_X, X]=hist(kvot,bins);
    
    hold on
    plot([1, X+.5*[diff(X),0]], 1- [0, cumsum(N_X)]/sum(N_X), '-', 'linewidth', 2)
    %plot(X, N_X/sum(N_X), '*');hold on
    data(:,2*i-1)=[0, X+.5*[diff(X),0]];
    data(:,2*i)=1- [0, cumsum(N_X)]/sum(N_X);
    
    
end

l=legend(leg_str);set(l, 'fontsize', 20)
xlabel('$\lambda_1/\lambda_2$','interpreter', 'latex')
ylabel('$1-\mathrm{CDF}$','interpreter', 'latex')

set(gca, 'fontsize', 20, 'yscale', 'log', 'xscale', 'lin')

% % save('egenvardeskvot.tsv', 'data', '-ascii')



%% Plottar ALLT, asymmetri
clc;clearvars;clf


b=dir('kvot*.mat');
leg_str=cell(length(b),1);
L=length(b);
bins=100;%logspace(0,2,50);

data=zeros(bins+1, 2*L);

for i=1:L
    leg_str{i}=b(i).name;
    load(leg_str{i}, 'egen')
    kvot=(egen(:,2)-egen(:,1)).^2./(egen(:,2)+egen(:,1)).^2;%assymetri
    [N_X, X]=hist(kvot,bins);
    hold on
    
    plot([0, X+.5*[diff(X),0]], 1- [0, cumsum(N_X)]/sum(N_X), '-', 'linewidth', 2)
    %plot([1, X+.5*[diff(X),0]], 1- [0, cumsum(N_X)]/sum(N_X), '-', 'linewidth', 2)
    %plot(X, N_X/sum(N_X), '*');hold on
    data(:,2*i-1)=[0, X+.5*[diff(X),0]];
    data(:,2*i)=1- [0, cumsum(N_X)]/sum(N_X);
    
    
end

l=legend(leg_str);set(l, 'fontsize', 20)
xlabel('$(\lambda_1-\lambda_2)^2/(\lambda_1+\lambda_2)^2$','interpreter', 'latex')
ylabel('$1-\mathrm{CDF}$','interpreter', 'latex')

set(gca, 'fontsize', 20, 'yscale', 'lin', 'xscale', 'lin')

% % save('asymmetri.tsv', 'data', '-ascii')








