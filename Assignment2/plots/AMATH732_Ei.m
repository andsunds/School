%%
clc;clf;clear

x=[5,10];
N=60;
Sn=zeros(length(x),N);
A=[2.5, 10.5];
data=zeros(N,5);

n=1:N;
data(:,1)=n';
symb='oxd*';
for j=1:2;
    a=A(j);
for i=n
    Sn(:,i)=S_Ei(i,x,a-1).';
end
data(:,[2*j 2*j+1])=log(abs(Sn).');
plot(n,abs(Sn(1,:)), ['-',symb(2*j-1)], 'markersize', 15)
hold on
plot(n,abs(Sn(2,:)), ['-',symb(2*j)], 'markersize', 15)

end

xlabel('$n$', 'interpreter','latex', 'fontsize',30)
ylabel('$|S_n|$', 'interpreter','latex', 'fontsize',30)

l=legend('$x=5,\phantom{0} a=2.5$','$x=10, a=2.5$','$x=5,\phantom{0} a=10.5$','$x=10, a=10.5$', 'location', 'NorthWest');
set(l,'interpreter','latex', 'fontsize',30)

set(gca,'yscale','log', 'fontsize',30)

axis([0,N,1e-5,1e20])