%%
clc;clf;clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

for fil=1:2 %för att plotta datan för båda celltyperna

data =load(filnamn{fil});


C = separera(data);
n=length(C);%antal partiklar

storl=zeros(n,1);%stoleken på partiklarna
medelsteg=zeros(n,1);%medelstegen på en prtikel
steg_osakerhet=zeros(n,1);%avvikelen från medelsteglängden
rorlighet=zeros(n,1);%hur stort område vandrar partikeln runt i

for i=1:n
    storl(i)=mean(C{i}(:,4));
    
    steglangd=sqrt(sum(diff(C{i}(:,2:3)).^2,2));%alla steglängder
    medelsteg(i)=mean(steglangd);
    steg_osakerhet(i)=std(steglangd)/sqrt(length(steglangd));
    
    rorlighet(i)=sqrt(sum(var(C{i}(:,2:3))));
end

figure(2*fil-1)
plot(storl, rorlighet, '.r')
title(filnamn{fil})
xlabel('Storlek', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('R\"o{}rlighet', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15);


figure(2*fil)
errorbar(storl, medelsteg, steg_osakerhet, '.')
title(filnamn{fil})
xlabel('Storlek', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('Stegl\"a{}ngd', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15);
axis([0,20,0,0.05])
end
