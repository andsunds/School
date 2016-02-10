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



%% varians som fkn av tid (Radius of duration)
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
STD=zeros(1000,1);
%STD_n=zeros(1000,1);
%loop över alla partiklar med tidsdata upp till 10s
tic
for i=find(cellfun('length',C)==1000).'; 
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    
    tmp=zeros(1000,2);%initialisering
    for j=1:1000 %loop över att tidpunkter
        tmp(j, :)=std(TN(1:j,:), 0, 1);
    end
    
    %bygger "medelvärde"
    STD=STD+tmp/(koef(1)*I(i).^koef(2));%normerat medelvärde
    %STD_n=STD_n+tmp(:,2)/(koef(1)*I(i).^koef(2));%normerat medelvärde
end
toc


%anpassar exponentialsamband, funkar inte
t=C{i}(2:end,1);
c=[ones(size(t)), log(t)]\log(STD(2:end,:).^2);

x=logspace( -4, log10(t(end)) ).';
y=exp(c(1,:)).*(x.^(c(2,:)));


%plottar data
subplot(1,2,fil)
plot(C{i}(:,1), STD.^2), hold on
plot(x,y)

str1=sprintf('%.1d t^{%1.2f}', exp(c(1,1)), c(2,1));
str2=sprintf('%.1d t^{%1.2f}', exp(c(1,2)), c(2,2));


legend('T', 'N', str1, str2, 'location', 'NorthWest')
title(filnamn{fil}(1:end-4))
xlabel('Tid/[s]', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('R\"o{}rlighet', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15)%,'XScale','log','YScale','log');
pause(.1)
end

%% S(dt)=(1/T) sum((f(t)-f(t+dt)).^2) over all t
clc;
figure(3);clf;pause(.1)

Dt=(1:300).';

for fil=1:2;
data =load(filnamn{fil});
C = separera(data);

%kortare namn
I=intensitet{fil};
lambda_t=rorlighet_t{fil};
lambda_n=rorlighet_n{fil};

%Detta är just nu samma normering som för rörligheten...
%Sammanvägd koef för T och N ger rörligheter som går att jämföra
lambda=sqrt(rorlighet_n{fil}.^2+rorlighet_t{fil}.^2);
%beräkna koefficienter för lutning i loglog
koef=[ones(size(I)), log(I)]\log(lambda);
koef(1)=exp(koef(1));%konverterar till potenssamband


S=zeros(length(Dt),2);
tic
for i=find(cellfun('length',C)==1000).'; 
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    
    tmp=zeros(length(Dt),2);
    for dt=Dt.'
        tmp(dt, :)=mean(diff(TN(1:dt:end,:), 1,1).^2, 1);
    end
%     f_t=@(dt) mean(diff(TN(1:dt:end,1), 1,1).^2, 1);
%     f_n=@(dt) mean(diff(TN(1:dt:end,2), 1,1).^2, 1);
%     tmp_t=arrayfun(f_t, Dt.');
%     tmp_n=arrayfun(f_n, Dt.');
%     tmp=[tmp_t, tmp_n];

    % Detta är samma normering som för rörligheten, 
    % kanske skulle man hitta på något annat.
    S=S+tmp/(koef(1)*I(i).^koef(2));
end
toc

%anpassar exponentialsamband
c=[ones(size(Dt)), log(Dt)]\log(S);

x=logspace(-1, log10(Dt(end))+0.05 ).';
y=exp(c(1,:)).*(x.^(c(2,:)));

%plottar data
subplot(1,2,fil)
plot(Dt,S), hold on
plot(x,y)

str1=sprintf('%.1d dt^{%1.2f}', exp(c(1,1)), c(2,1));
str2=sprintf('%.1d dt^{%1.2f}', exp(c(1,2)), c(2,2));

legend('T', 'N', str1, str2, 'location', 'NorthWest')
title(filnamn{fil}(1:end-4))
xlabel('Tid/[s]', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('R\"o{}rlighet', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15)%,'XScale','log','YScale','log');
pause(.1)
end


