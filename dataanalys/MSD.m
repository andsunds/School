%MSD

%MSD for the small particles, quite linear
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