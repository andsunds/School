%CTRW
clf
%F?rs?k till att simulera CTRW
n=1;% Antal partiklar
N=900000;%Anpassas efter hur alpha v?ljs
N_trials=1e5;

%Generera v?ntetider, lite hokus pokus och sopande under mattan
F_i=@(U,eps,alpha,A)(eps^(-alpha)-alpha*U/A).^(-1/alpha); %Invers CDF, CDF=tau^-(1+alpha)
alpha=0.65; %Anpassningsbar parameter
eps=1e-7; %Minsta m?jliga v?ntetid
A=alpha*eps^alpha; %Normeringskonstant

kvot=zeros(N_trials,1);
egen=zeros(N_trials,2);

tic

for k=1:N_trials
    
U=rand(N,n);
tau=F_i(U,eps,alpha,A); %V?ntetider
t=cumsum(tau,1); %tid

X_CTRW = randn(N,2,n);%Normalf?rdelade steg

XY=cumsum(X_CTRW,1); %x,y separerade
%R=sqrt(XY(:,1).^2+XY(:,2).^2);
%subplot(2,1,1)
%plot(t,XY)
%xlabel('tid (s)')
%ylabel('F?rflyttning')
%title('XY')
%subplot(2,1,2)
%plot(t,R)
%xlabel('tid (s)')
%ylabel('F?rflyttning')
%title('R')


%Asymmetri

%Sampla med jämna tidssteg
N_sampl=1000;
t_even=linspace(0,1e1,N_sampl);
t_begin=1;
XY_even=zeros(N_sampl,2);
for i=1:N_sampl
    for j=t_begin:N
        if (t_even(i)<=t(j))
            XY_even(i,:)=XY(j,:);
            j_begin=j;
            break
        end
    end
end

%subplot(2,1,1)
%hold on
%plot(t_even,XY_even(:,1))
%hold off

%Asphericity
    X=XY_even(:,1);
    Y=XY_even(:,2);
    gyr = [sum(X.^2)/N_sampl-(sum(X)/N_sampl)^2, sum(X.*Y)/N_sampl-(sum(Y)/N_sampl)*(sum(X)/N_sampl); 
           sum(Y.*X)/N_sampl-(sum(Y)/N_sampl)*(sum(X)/N_sampl), sum(Y.^2)/N_sampl-(sum(Y)/N_sampl).^2];
       
    [~,D] = eig(gyr);
    %D(1)
    %D(4)
    kvot(k)=( (D(4)-D(1))/(D(4)+D(1)) )^2;
    egen(k,:)=[D(1), D(4)];

end

%kvot
toc

%hist(kvot)

bins=200;%logspace(0,2);
[N_X, X]=hist(kvot,bins);
plot([0, X+.5*[diff(X),0]], 1- [0, cumsum(N_X)]/sum(N_X), '-')
set(gca, 'fontsize', 20, 'yscale', 'lin', 'xscale', 'lin')

%%
%Plotta ~58k simulerade CTRW

egen_CTRW=load('egenCTRW58k.tsv');
kvot_CTRW=(abs(egen_CTRW(:,1)-egen_CTRW(:,2)))./(egen_CTRW(:,1)+egen_CTRW(:,2));

bins=200;%logspace(0,2);
[N_X, X]=hist(kvot_CTRW,bins);
plot([0, X+.5*[diff(X),0]], 1- [0, cumsum(N_X)]/sum(N_X), '-')
set(gca, 'fontsize', 20, 'yscale', 'lin', 'xscale', 'lin')