% Hur många partiklar av olika storlek
figure(2), clf, figure(3), clf, figure(1), clf
clear

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';
filstr={'energydepleted' 'logphase'};%For usage in titles

for k=1:2 %For both files
fil=k;
subplot(1,2,k)
data =load(filnamn{fil});

C = separera(data);
n=length(C);%Number of particles

intensitet=zeros(1,n); %Intensity of the particles
for i=1:n
    intensitet(i)=C{i}(1,4);
end



%Intensity distribution
figure(1)
subplot(1,2,k)

bins=20;
[bin_counts bin_center]=hist(intensitet,bins); %Number of particles per box
hist(intensitet,bins)
str_k=strcat('Intensitetsfördelning, ',filstr{k});
title(str_k)

%To compare the size of the smallest to the largest particle
intensitetskvot=max(intensitet)/min(intensitet);
areakvot(k)=(intensitetskvot)^(2/3); %Area ratio between largest and smallest particle



%Sort particles by increasing intensity into boxes with particles of
%similar intensity in the same box
[sorterade_intensiteter index_innan]=sort(intensitet); %To get the indeces for the particles in the histogram boxes

index_particles_sizesorted=cell(1,bins); %Indeces for particles belonging to the same intensitybox
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
  


% MSD for the first boxes (lower intensity, more datapoints)
figure(2)
subplot(1,2,k)


j_max=7; %Number of bins used in calculation
n_max=916; %Number of data points (time-position), 916=lowest index
MSD_tot=0;

for j=1:j_max
    MSD=zeros(1,n_max)';
    for i=index_particles_sizesorted{j};
        XY = C{i}(:,2:3);
        Intensitet_i=C{i}(1,4);
        MSD_i=sum(XY.^2,2)*(Intensitet_i)^(1/3); %To remove intensity dependency
        MSD=MSD+MSD_i(1:n_max);
    end
    MSD=MSD/bin_counts(j); %MSD mean between similar intensities
    MSD_tot=MSD_tot+MSD*bin_counts(j); %Weight by number of particles
    hold on
    plot(C{1}(1:n_max,1),MSD)
end
hold off
center_etikett=num2str(bin_center,3);
title(strcat('MSD, boxar, ',filstr{k}))
legend(center_etikett(1:5), center_etikett(6+6:(12+4)), center_etikett((16+6):(22+4)), center_etikett((26+6):(32+4)),center_etikett((36+6):(42+4)), center_etikett(46+6:(52+4)), center_etikett(56+6:(62+4)),'Location','Best')
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


%Finds exponential fit y=a*exp(b*x)

show=901; %Number of timesteps
DT=(1:n_max)';
Dt=(DT(1:show)-1)*1e-2; %real time
start=2;
c=[ones(show+1-start,1), log(Dt(start:end))]\log(MSD_tot(start:show,:)); %To find constants a,b

c_both(k,1:2)=c; %Store constants for both plots

x=logspace(-2, log10(Dt(end)) ).';
y=repmat(exp(c(1,:)), length(x),1).* bsxfun(@power, x, c(2,:));

%plots fitted graph
hold on
plot(Dt,MSD_tot(1:show,:)), hold on
plot(x,y)
title(strcat('MSD, viktad summa,',filstr{k}))
xlabel('tid')
ylabel('MSD')

str1=sprintf('%.1d dt^{%1.2f}', exp(c(1,1)), c(2,1)); %FIELD WIDTH.PRECISION SUBTYPE, subtypes: d=signed integer f=fixedpoint float
legend(str1, 'location', 'North')
end
