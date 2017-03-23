%%
clc;clear;
figure(1)
clf
figure(2)
clf


R=@(n1,n2,n3, D) ...
    ( (n2*(n1-n3)*cos(D*n2)).^2 + ((n2^2-n1*n3)*sin(D*n2)).^2 )...
    ./( (n2*(n1+n3)*cos(D*n2)).^2 + ((n2^2 + n1*n3)*sin(D*n2)).^2 );
   

T=@(n1,n2,n3, D) ...
    ( 4*n1.^2*n2.^2 )...
    ./( (n2*(n1+n3)*cos(D*n2)).^2 + ((n2^2 + n1*n3)*sin(D*n2)).^2 );



D0=linspace(0,1,1000); %Use D0=d/lambda, so in R and T use 2*pi*D.

n1=1;
n2=[1, 1.2, 3];
n3=1.5;
linetype = {'-k', '--k', '-.k'};

for i=1:3
figure(1)
hold on
plot(D0, R(n1,n2(i),n3, 2*pi*D0), linetype{i})
figure(2)
hold on
plot(D0, T(n1,n2(i),n3, 2*pi*D0), linetype{i})
     %D0, R(n1,n2,n3, 2*pi*D0)+T(n1,n2,n3, 2*pi*D0))
end


l=legend(sprintf('$n_2=%1.1f$',n2(1)),...
    sprintf('$n_2=%1.1f$',n2(2)),sprintf('$n_2=%1.1f$',n2(3)));
set(l, 'interpreter', 'latex')
xlabel('$d/\lambda_0$','interpreter', 'latex')
ylabel('$T(d)$','interpreter', 'latex')
set(gca, 'fontsize', 12, 'ylim', [0,1])


figure(1)
l=legend(sprintf('$n_2=%1.1f$',n2(1)),...
    sprintf('$n_2=%1.1f$',n2(2)),sprintf('$n_2=%1.1f$',n2(3)));
set(l, 'interpreter', 'latex')
xlabel('$d/\lambda_0$','interpreter', 'latex')
ylabel('$R(d)$','interpreter', 'latex')
set(gca, 'fontsize', 12, 'ylim', [0,1])


























