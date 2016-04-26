%% Steglängd och rörlighet(std) som fkn av storlek
clc;
figure(1);clf;pause(.1)
clearvars

load('filnamn.mat')

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



% % save('kompleterande_data.mat', 'intensitet', 'medelsteg', 'steg_std', 'std_t', 'std_n', 'sigma_brus', '-mat')

%% Spara för plottning
clc;clearvars

load('filnamn.mat')
load(kompl)

n=[length(intensitet{1}),length(intensitet{2})];
plot_data=NaN(max(n), 6);

for fil=1:2
plot_data(1:n(fil),(3*fil-2):(3*fil))=[intensitet{fil}, sqrt(std_n{fil}.^2+std_t{fil}.^2), medelsteg{fil}];
end

save('storleksberoende.tsv', 'plot_data', '-ascii')





