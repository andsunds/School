%%
clc;clf;clear

x=[5,10];
N=60;
Sn=zeros(length(x),N);

n=1:N;
for a=[2.5, 10.5];
for i=n
    Sn(:,i)=S_Ei(i,x,a-1).';
end

plot(n,abs(Sn), '-o')
hold on
end

xlabel('$n$', 'interpreter','latex', 'fontsize',30)
ylabel('$|S_n|$', 'interpreter','latex', 'fontsize',30)

l=legend('$x=5, a=2.5$','$x=10, a=2.5$','$x=5, a=10.5$','$x=10, a=10.5$', 'location', 'NorthWest');
set(l,'interpreter','latex', 'fontsize',30)

set(gca,'yscale','log')