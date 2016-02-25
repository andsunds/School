% Hur många partiklar av olika storlek
figure(2), clf, figure(3), clf, figure(1), clf
clear

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';
filstr={'energydepleted' 'logphase'};%For usage in titles

for k=1:2%
fil=k;
subplot(1,2,k)
data =load(filnamn{fil});

C = separera(data);
n=length(C);%antal partiklar
%index_particle_0_128=[]; %Index för partiklar upp till intensiteten 1.28
%index_particle_129_237=[];
%index_largeparticle=[]; %Index för partiklar med intensitet mellan 8.9-21.8

intensitet=zeros(1,n);
for i=1:n
    intensitet(i)=C{i}(1,4);
    %if (intensitet(i)<1.28)
    %    index_particle_0_128=[index_particle_0_128 i];
    %elseif (8.9<intensitet(i))
    %    index_largeparticle=[index_largeparticle i];
    %end
end
%Hela denna for-loop kan ersättas med:
%index_smallparticle=find(intensitet(i)<1.28));
%index_largeparticle=find(8.9<intensitet(i));

figure(1)

bins=20;
[bin_counts bin_center]=hist(intensitet,bins); %Antal per låda
hist(intensitet,bins)
str_k=strcat('Intensitetsfördelning, ',filstr{k});
title(str_k)

intensitetskvot=max(intensitet)/min(intensitet);
areakvot(k)=(intensitetskvot)^(2/3); %En faktor ~20 skiljer i areastorleken mellan de största och minsta partiklarna i logphase

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
  

% MSD for the smaller particles

figure(2)
subplot(1,2,k)


j_max=5; %Number of bins used in calculation
n_max=916; %Number of data points
MSD_tot=0;

for j=1:j_max %length(index_particles_sizesorted)
    MSD=zeros(1,n_max)'; %916=lowest index
    for i=index_particles_sizesorted{j};
        XY = C{i}(:,2:3);
        Intensitet_i=C{i}(1,4)
        MSD_i=sum(XY.^2,2)*(Intensitet_i)^(1/3); %To remove intensity dependency
        MSD=MSD+MSD_i(1:n_max);
    end
    MSD=MSD/bin_counts(j);%;
    MSD_tot=MSD_tot+MSD*bin_counts(j); %Weight by number of particles
    hold on
    plot(C{1}(1:n_max,1),MSD)
end
hold off
center_etikett=num2str(bin_center,3);
title(strcat('MSD, boxar, ',filstr{k}))
legend(center_etikett(1:5), center_etikett(6+6:(12+4)), center_etikett((16+6):(22+4)), center_etikett((26+6):(32+4)),center_etikett((36+6):(42+4)),'Location','Best')
xlabel('tid')
ylabel('MSD')

figure(3)
subplot(1,2,k)

MSD_tot=MSD_tot/sum(bin_counts(1:j_max)); %Divide by total number of particles

%Plottar i annan tidsskala än nedan
%plot(C{1}(1:n_max,1),(MSD_tot)) 
%title(strcat('MSD, viktad summa,',filstr{k}))
%xlabel('tid')
%ylabel('MSD')


%anpassar exponentialsamband

show=901; %hur många tidssteg ska undersökas
DT=(1:n_max)';
%MD_tot=(MSD_tot); %Mean displacement
Dt=(DT(1:show)-1)*1e-2;%verklig tid
start=2;
c=[ones(show+1-start,1), log(Dt(start:end))]\log(MSD_tot(start:show,:)); %

c_both(k,1:2)=c;

x=logspace(-2, log10(Dt(end)) ).';
y=repmat(exp(c(1,:)), length(x),1).* bsxfun(@power, x, c(2,:));

hold on
%plottar anpassad kurva
plot(Dt,MSD_tot(1:show,:)), hold on
plot(x,y)
title(strcat('MSD, viktad summa,',filstr{k}))
xlabel('tid')
ylabel('MSD')

str1=sprintf('%.1d dt^{%1.2f}', exp(c(1,1)), c(2,1))
legend(str1, 'location', 'North')
end
