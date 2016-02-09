%% Steglängd och rörlighet(std) som fkn av storlek
clc;
figure(1);clf;
clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

%cellinitialisering
intensitet=cell(1,2);
medelsteg=cell(1,2);
steg_osakerhet=cell(1,2);
rorlighet_t=cell(1,2);
rorlighet_n=cell(1,2);


for fil=1:2 %för att undersöka datan för båda celltyperna
%laddar in data
data =load(filnamn{fil});
C = separera(data);
n=length(C);%antal partiklar

%Initialisering av enskilda celler
intensitet{fil}=zeros(n,1);%stoleken på partiklarna
medelsteg{fil}=zeros(n,1);%medelstegen på en prtikel
steg_osakerhet{fil}=zeros(n,1);%avvikelen från medelsteglängden
%L=zeros(n,1);U=zeros(n,1);
rorlighet_t{fil}=zeros(n,1);%hur stort område vandrar partikeln runt i
rorlighet_n{fil}=zeros(n,1);%hur stort område vandrar partikeln runt i

for i=1:n
    intensitet{fil}(i)=mean(C{i}(:,4));
    
    TN=koordinatbyte(C{i}(:,2:3));%positionsdata i TN-koordinater
    
    steglangd=sqrt(sum(diff(TN).^2, 2));%alla steglängder
    medelsteg{fil}(i)=mean(steglangd);
    steg_osakerhet{fil}(i)=std(steglangd);%/sqrt(length(steglangd));
    
    rorlighet_t{fil}(i)=sqrt(var(TN(:,1)));
    rorlighet_n{fil}(i)=sqrt(var(TN(:,2)));
end
%plottar rörrilghet
subplot(2,2,2*fil-1)
plot(intensitet{fil}, rorlighet_t{fil}, 'or'), hold on
plot(intensitet{fil}, rorlighet_n{fil}, 'xg')
title(filnamn{fil}(1:end-4))
xlabel('Intensitet', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('R\"o{}rlighet', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15)%,'XScale','log','YScale','log');


subplot(2,2,2*fil)
errorbar(intensitet{fil}, medelsteg{fil}, steg_osakerhet{fil}, '.')
title(filnamn{fil}(1:end-4))
xlabel('Intensitet', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('Stegl\"a{}ngd', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15)%,'XScale','log','YScale','log');
axis([0,20,0,8e-9])

%pause(0.1)%för att rita ut efter hand


% %Hittar potenssamband
% fprintf('\n%s:',filnamn{fil}(1:end-4))
% fprintf('\n\nPotens längs T:')
% fit(intensitet{fil},rorlighet_t{fil},'power1') %Powerfit r�rlighet
% fprintf('\n\nPotens längs N:')
% fit(intensitet{fil},rorlighet_n{fil},'power1') %Powerfit r�rlighet
% fprintf('\n\n')

end



%% varians som fkn av tid
clc;
figure(2);clf;pause(.1)

for fil=1:2;
data =load(filnamn{fil});
C = separera(data);

%kortare namn
I=intensitet{fil};
lambda_t=rorlighet_t{fil};
lambda_n=rorlighet_n{fil};

%beräkna koefficienter för lutning i loglog
% %Olika koef för T och N leder till samma rörlighet i slutet
% koef_t=[ones(size(I)), log(I)]\log(lambda_t);
% koef_t(1)=exp(koef_t(1));%konverterar till potenssamband
% koef_n=[ones(size(I)), log(I)]\log(lambda_n);
% koef_n(1)=exp(koef_n(1));%konverterar till potenssamband

%Sammanvägd koef för T och N ger rörligheter som går att jämföra
lambda=sqrt(rorlighet_n{fil}.^2+rorlighet_t{fil}.^2);
%beräkna koefficienter för lutning i loglog
koef=[ones(size(I)), log(I)]\log(lambda);
koef(1)=exp(koef(1));%konverterar till potenssamband


% figure(1) %Ritar in anpassning
% subplot(2,2,2*fil-1)
% x=logspace(-4, 1, 1000)*2;
% plot(x,koef_t(1)*x.^koef_t(2),'-k')
% plot(x,koef_n(1)*x.^koef_n(2),'--k')
% figure(2)

%initialisering
STD_t=zeros(1000,1);
STD_n=zeros(1000,1);
%loop över alla partiklar med tidsdata upp till 10ms
for i=find(cellfun('length',C)==1000).'; 
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    
    tmp=zeros(1000,2);%initialisering
    for j=1:1000 %loop över att tidpunkter
        tmp(j, :)=std(TN(1:j,:), 0, 1);
    end
    %bygger "medelvärde"
    %STD_t=STD_t+tmp(:,1)/(koef_t(1)*I(i).^koef_t(2));%normerat medelvärde
    %STD_n=STD_n+tmp(:,2)/(koef_n(1)*I(i).^koef_n(2));%normerat medelvärde
    %bygger "medelvärde"
    STD_t=STD_t+tmp(:,1)/(koef(1)*I(i).^koef(2));%normerat medelvärde
    STD_n=STD_n+tmp(:,2)/(koef(1)*I(i).^koef(2));%normerat medelvärde
end

%plottar data
subplot(1,2,fil)
plot(C{i}(:,1), STD_t), hold on
plot(C{i}(:,1), STD_n)
legend('Rörlighet T', 'Rörlighet N', 'location', 'SouthEast')
title(filnamn{fil}(1:end-4))
xlabel('Tid/[s]', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('R\"o{}rlighet', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15,'XScale','log','YScale','log');
pause(.1)
end


%% Hur många partiklar av olika storlek
clf

fil=2;
data =load(filnamn{fil});

C = separera(data);
n=length(C);%antal partiklar
index_smallparticle=[]; %Index för partiklar upp till intensiteten 1.28
index_largeparticle=[]; %Index för partiklar med intensitet mellan 8.9-21.8

intensitet=zeros(1,n);
for i=1:n
    intensitet(i)=C{i}(1,4);
    if (intensitet(i)<1.28)
        index_smallparticle=[index_smallparticle i];
    elseif (8.9<intensitet(i))
        index_largeparticle=[index_largeparticle i];
    end
end%Hela denna for-loop kan ersättas med:
%index_smallparticle=find(intensitet(i)<1.28));
%index_largeparticle=find(8.9<intensitet(i));

bins=40;
hist(intensitet,bins)

intensitetskvot=max(intensitet)/min(intensitet);
areakvot=(intensitetskvot)^(2/3) %

%En faktor ~20 skiljer i areastorleken mellan de största och minsta partiklarna
