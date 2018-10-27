%% Absorption
clc;clear all;
figure(3), clf

load('Data/halogen_backgound_absorption.mat')

col=cell(1,2);
col{1}='pink';
col{2}='green';

plotcolor=cell(1,2);
plotcolor{1}='-.m';
plotcolor{2}='-g';


for j=1:2

a=dir(['Data/halogen_', col{j} ,'stuff_kompl*']);
b=dir(['Data/halogen_', col{j} ,'stuff_0*']);

La=length(a);%olika långa
Lb=length(b);%olika långa

Xa=zeros(211,La); Xb=zeros(701,Lb);
Ya=zeros(211,La); Yb=zeros(701,Lb);


for i=1:La
    data_a=load(a(i).name);
    Xa(:,i)=data_a(:,2);
    Ya(:,i)=data_a(:,1);
end

for i=1:Lb
    data_b=load(b(i).name);
    Xb(:,i)=data_b(:,2);
    Yb(:,i)=data_b(:,1);
end

la=mean(Xa,2); lb=mean(Xb,2);
Ia=mean(Ya,2); Ib=mean(Yb,2);



L=[la(1:200); (la(201)+lb(1:10))/2; lb(11:end)];
I=[Ia(1:200); (Ia(201)+Ib(1:10))/2; Ib(11:end)];
korr=exp(polyval(P_korr, L));



% figure(1), clf
% plot(L,I.*korr)
% hold on
% plot(L,I_background.*korr)
% axis([L_min, 1100, 1e-2, max(I_background(L>L_min).*korr(L>L_min))*1.1])
% legend(sprintf('%s stuff', col), 'background', 'location', 'best')
% 
% set(gca, 'fontsize', 15,'yscale', 'log')
% 
% ylabel('Intensitet /[a.u.]', 'interpreter', 'latex', 'fontsize', 20)
% xlabel('$\lambda$/[nm]', 'interpreter', 'latex', 'fontsize', 20)



figure(3)
A=(I./I_background);
plot(L,A, plotcolor{j})
hold on
end

axis([300, 1100, .2, 1.15])

set(gca, 'fontsize', 15, 'yscale', 'log')

l=legend('Rhodamin B', 'Kumarin 307', 'location', 'best');
set(l, 'interpreter', 'latex', 'fontsize', 15)

ylabel('$I_\mathrm{abs.}/I_\mathrm{bakgrund}$', 'interpreter', 'latex', 'fontsize', 20)
xlabel('$\lambda$/[nm]', 'interpreter', 'latex', 'fontsize', 20)

