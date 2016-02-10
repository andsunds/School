% Hur många partiklar av olika storlek
clf, clear

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=2;
data =load(filnamn{fil});

C = separera(data);
n=length(C);%antal partiklar
index_particle_0_128=[]; %Index för partiklar upp till intensiteten 1.28
index_particle_129_237=[];
index_largeparticle=[]; %Index för partiklar med intensitet mellan 8.9-21.8

intensitet=zeros(1,n);
for i=1:n
    intensitet(i)=C{i}(1,4);
    if (intensitet(i)<1.28)
        index_particle_0_128=[index_particle_0_128 i];
    elseif (8.9<intensitet(i))
        index_largeparticle=[index_largeparticle i];
    end
end
%Hela denna for-loop kan ersättas med:
%index_smallparticle=find(intensitet(i)<1.28));
%index_largeparticle=find(8.9<intensitet(i));

bins=20;
[bin_counts bin_center]=hist(intensitet,bins); %Antal per låda
hist(intensitet,bins)

intensitetskvot=max(intensitet)/min(intensitet);
areakvot=(intensitetskvot)^(2/3); %En faktor ~20 skiljer i areastorleken mellan de största och minsta partiklarna

[sorterade_intensiteter index_innan]=sort(intensitet); %För att plocka ut index från histogramboxarna

index_particles_sizesorted=cell(1,bins);
index_particles_sizesorted{1}=index_innan(1:bin_counts(1));
total=bin_counts(1);
for i=2:bins
    if bin_counts==0
        index_particles_sizesorted{i}=0;
    else
        index_particles_sizesorted{i}=index_innan((total+1):((total+1)+bin_counts(i)-1));
    end
    total=total+bin_counts(i);
end
    


%% MSD for the small particles
clf;


fil=2;
data = separera(load(filnamn{fil}));

for j=1:5 %length(index_particles_sizesorted)
    MSD=zeros(1,916)'; %916=l?gsta index
    for i=index_particles_sizesorted{j};
        XY = data{i}(:,2:3);
        MSD_i=sum(XY.^2,2);
        MSD=MSD+MSD_i(1:916);
    end
    MSD=MSD/length(index_particles_sizesorted{j})
    hold on
    plot(data{1}(1:916,1),MSD)
end

center_etikett=num2str(bin_center,3)
title('MSD, olika intensitet')
legend(center_etikett(1:5), center_etikett(6+6:(12+4)), center_etikett((16+6):(22+4)), center_etikett((26+6):(32+4)),center_etikett((36+6):(42+4)),'Location','Best')
xlabel('tid')
ylabel('MSD')