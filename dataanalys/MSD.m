
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

%% MSD for the small particles, quite linear
clf;

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=2;
data = separera(load(filnamn{fil}));

MSD=zeros(1,916)'; %916=l?gsta index

for i=index_smallparticle;
    XY = data{i}(:,2:3);
    MSD_i=sum(XY.^2,2);
    MSD=MSD+MSD_i(1:916);
end

MSD=MSD/length(index_smallparticle)

plot(data{1}(1:916,1),MSD)
title('MSD small particles')