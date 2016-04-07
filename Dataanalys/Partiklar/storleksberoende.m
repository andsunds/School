%% Steglängd och rörlighet(std) som fkn av storlek
clc;
figure(1);clf;pause(.1)
clearvars

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

%cellinitialisering
intensitet=cell(2,1);
medelsteg=cell(2,1);
steg_std=cell(2,1);
std_t=cell(2,1);
std_n=cell(2,1);
sigma_brus=zeros(2,1);


for fil=1:2 %för att undersöka datan för båda celltyperna
%laddar in data
data =load(filnamn{fil});
C = separera(data);
n=length(C);%antal partiklar

%Initialisering av enskilda celler
intensitet{fil}=zeros(n,1);%stoleken på partiklarna
medelsteg{fil}=zeros(n,1);%medelstegen på en prtikel
steg_std{fil}=zeros(n,1);%avvikelen från medelsteglängden

std_t{fil}=zeros(n,1);%hur stort område vandrar partikeln runt i
std_n{fil}=zeros(n,1);%hur stort område vandrar partikeln runt i

for i=1:n
    intensitet{fil}(i)=mean(C{i}(:,4));
    
    TN=koordinatbyte(C{i}(:,2:3));%positionsdata i TN-koordinater
    
    steglangd=(sum(diff(TN).^2, 2));%alla steglängder
    medelsteg{fil}(i)=sqrt(mean(steglangd));
    steg_std{fil}(i)=std(steglangd);
    
    std_t{fil}(i)=sqrt(var(TN(:,1)));
    std_n{fil}(i)=sqrt(var(TN(:,2)));
end

index=find(intensitet{fil}>12);
sigma_brus(fil) = 0.5*mean( medelsteg{fil}(index) );

%plottar rörrilghet
subplot(2,2,2*fil-1)
%plot(intensitet{fil}, rorlighet_t{fil}, 'or'), hold on
%plot(intensitet{fil}, rorlighet_n{fil}, 'xg')
plot(intensitet{fil}, sqrt(std_n{fil}.^2+std_t{fil}.^2), 'k*')
hold on
plot([0,25], sqrt(2)*sigma_brus(fil)*[1,1], 'r');hold off

title(filnamn{fil}(1:end-4))
xlabel('Intensitet', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('Standardavvikelse /[m]', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15,'XScale','lin','YScale','log');
axis([0,25,5e-10,10e-8])

subplot(2,2,2*fil)
%errorbar(intensitet{fil}, medelsteg{fil}, steg_osakerhet{fil}, '.')

plot(intensitet{fil}, medelsteg{fil}, '.'); hold on
plot([0,25], 2*sigma_brus(fil)*[1,1], 'r');hold off
title(filnamn{fil}(1:end-4))
xlabel('Intensitet', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('Stegl\"a{}ngd /[m]', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15, 'ylim', [0,4e-17],'XScale','lin','YScale','log');
axis([0,25,7e-10,1e-8])

%pause(0.1)%för att rita ut efter hand

end


% % save('kompleterande_data.mat', 'filnamn', 'intensitet', 'medelsteg', 'steg_std', 'std_t', 'std_n', 'sigma_brus', '-mat')


%% s(dt) = 1/'#particles' * sum( ((x(t)-x(0)).^2 ) over all particles
clc; clearvars
figure(2);clf;pause(.1)

load('kompleterande_data.mat',...
     'filnamn', 'intensitet', 'std_n', 'std_t', 'sigma_brus', '-mat')

for fil=1:2;
data =load(filnamn{fil});
C = separera(data);

%kortare namn
I=intensitet{fil};


%Sammanvägd koef för T och N ger rörligheter som går att jämföra
lambda=sqrt(std_n{fil}.^2+std_t{fil}.^2 - 2*sigma_brus(fil)^2 );
%beräkna koefficienter för lutning i loglog
koef=[ones(size(I)), log(I)]\log(lambda);
koef(1)=exp(koef(1));%konverterar till potenssamband




% figure(1) %Ritakoef = [1, -1/3];r in anpassning
% subplot(2,2,2*fil-1)
% x=logspace(-4, 1, 1000)*2;
% plot(x,koef_t(1)*x.^koef_t(2),'-k')
% plot(x,koef_n(1)*x.^koef_n(2),'--k')
% figure(2)

%initialisering
s=zeros(1000,2);
%STD_n=zeros(1000,1);
%loop över alla partiklar med tidsdata upp till 10s
tic
index=find(cellfun('length',C)==1000).'; 
for i=index;
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    %{
    tmp=zeros(1000,2);%initialisering
    for j=1:1000 %loop över at(koef(1)*I(i).^koef(2))t tidpunkter
        tmp(j, :)=var(TN(1:j,:), 0, 1);
    end
    %}
    tmp=(TN.^2);
    
    %bygger "medelvärde"
    s=s+tmp/(koef(1)*I(i).^koef(2));%normerat medelvärde
    
end
toc

s=sum(s,2)/length(index);


start=10;
stop=101;
t=C{i}(:,1);
c=[ones(stop+1-start,1), log(t(start:stop))]\log(s(start:stop,:));

x=logspace( -2, log10(t(end))).';
y=repmat(exp(c(1,:)), length(x),1).* bsxfun(@power, x, c(2,:));


%plottar data
subplot(1,2,fil)
plot(t, s), hold on
plot(x,y)

str1=sprintf('%.1d t^{%1.2f}', exp(c(1,1)), c(2,1));
%str2=sprintf('%.1d t^{%1.2f}', exp(c(1,2)), c(2,2));


%legend('T', 'N', str1, str2, 'location', 'NorthWest')
legend('Data', str1, 'location', 'NorthWest')
title(filnamn{fil}(1:end-4))
xlabel('Tid/[s]', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15,'XScale','log','YScale','log');
pause(.1)
end

%% S(dt)=(1/T) sum((f(t)-f(t+dt)).^2) over all t
clc; clearvars
figure(3);clf;
%figure(4);clf;
pause(.1)

load('kompleterande_data.mat',...
     'filnamn', 'intensitet', 'std_n', 'std_t', 'sigma_brus', '-mat')

N=1000;
DT=(1:N).';
addpath('../');%Lägger till så att create_indecis kan användas
INDECIES=create_indecis(N);
LENGHTS=N+1-DT;



for fil=1:2;
data =load(filnamn{fil});
C = separera(data);

%kortare namn
I=intensitet{fil};

%Detta är just nu samma normering som för rörligheten...
lambda=sqrt(std_n{fil}.^2+std_t{fil}.^2 -2*sigma_brus(fil).^2 );
%beräkna koefficienter för lutning i loglog
koef=[ones(size(I)), log(I)]\log(lambda);
koef(1)=exp(koef(1));%konverterar till potenssamband

%koef = [1, -1/3];

S=zeros(N,2);
tic
index=find(cellfun('length',C)==N).'; 
for i=index
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    
    tmp=zeros(length(DT),2);
    
    %Hög minnesåtgång
    A=repmat(TN(:,1), 1,N);
    A=triu(A.'-A);
    tmp(:,1)=sum(A(INDECIES).^2, 2)./LENGHTS;
    
    A=repmat(TN(:,2), 1,N);
    A=triu(A.'-A);
    tmp(:,2)=sum(A(INDECIES).^2, 2)./LENGHTS;
    
    
    % Detta är samma normering som för rörligheten, 
    % kanske skulle man hitta på något annat.
    S=S+(tmp)/(koef(1)*I(i).^koef(2));
end
toc

S=sum(S,2)/length(index);

show=101; %hur många tidssteg ska undersökas

%anpassar exponentialsamband
Dt=(DT-1)*1e-2;%verklig tid
start=10;
c=[ones(show+1-start,1), log(Dt(start:show))]\log(S(start:show,:));

x=logspace(-2, log10(Dt(end)) ).';
y=repmat(exp(c(1,:)), length(x),1).* bsxfun(@power, x, c(2,:));


% plottar data
figure(3);
subplot(1,2,fil)
plot(Dt,S), hold on
plot(x,y)

%uncertainty=1./sqrt(N+1-DT(1:show));
%plot(Dt, S(1:show,:).*(1+repmat(uncertainty,1,2)), '--k')
%q=plot(Dt, S(1:show,:).*(1-repmat(uncertainty,1,2)), '--k')
%legend(q, 'Osäkerhetsgränser')

%axis([0, Dt(end), 0, 6e-6]);

str1=sprintf('%.1d dt^{%1.2f}', exp(c(1,1)), c(2,1));
%str2=sprintf('%.1d dt^{%1.2f}', exp(c(1,2)), c(2,2));

%legend('T', 'N', str1,str2, 'location', 'NorthWest')%Om man seprarerar T och N
legend('T', str1, 'location', 'NorthWest')
title(filnamn{fil}(1:end-4))
xlabel('Tidssteg dt/[s]', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('S(dt)', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15,'XScale','log','YScale','log');
pause(.1)


end






