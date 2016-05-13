%fBm

clc;clf; clearvars
hold on
N_steps=1000;
N_trials=1000;
H=0.2;

%for brus=[0, 6];
%brus=0;
tic

%kvot=zeros(N_trials,1);
egen=zeros(N_trials,2);
for i=1:N_trials
    X=wfbm(H,N_steps); %+ brus*randn(N_steps,1);
    Y=wfbm(H,N_steps); %+ brus*randn(N_steps,1);
    X=X-mean(X);Y=Y-mean(Y);
    gyr = [sum(X.^2), sum(X.*Y); 
           sum(Y.*X), sum(Y.^2)]/N_steps;
    
    %gyr = [sum(X.^2)-sum(X).^2/N_steps, sum(X.*Y)-sum(X).*sum(Y)/N_steps; 
    %       sum(Y.*X)-sum(X).*sum(Y)/N_steps, sum(Y.^2)-sum(Y).^2/N_steps]/N_steps;
    [~,D] = eig(gyr);

    %kvot(i)=D(4)/D(1);%^((randn>0)*2 -1 );
    egen(i,:)=[D(1), D(4)];
end
%sparar sim

% % save(sprintf('egen_brus_%2.1f.mat', brus), 'egen', '-mat')

toc


bins=50;%logspace(0,2);
A=(egen(:,2)-egen(:,1)).^2./(egen(:,2)+egen(:,1)).^2;%assymetri
[N_X, X]=hist(A,bins);
plot([0, X+.5*[diff(X),0]], 1- [0, cumsum(N_X)]/sum(N_X), '-')
set(gca, 'fontsize', 20, 'yscale', 'lin', 'xscale', 'lin')
%end