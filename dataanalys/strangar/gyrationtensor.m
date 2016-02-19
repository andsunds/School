% Gyration tensor
clear all; clf;clc;
filnamn=cell(1,4);
filnamn{1}='confined_270304-6-28-min.mat';
filnamn{2}='confined_280204-2-32min.mat';
filnamn{3}='nonconfined_180304-1-5min.mat';
filnamn{4}='nonconfined_250104-1-167min.mat';

A = importdata(filnamn{1});
nbrB = size(A,3);%Antal bilder

XY=cell(nbrB,1); % Cell med positioner
G = cell(nbrB,1); % Gyration tensor

for i=1:nbrB
    [rad,kol] = find(255==A(:,:,i)); % Hitta index för position
    XY{i}=[kol,rad];
end

CM = zeros(2,1,nbrB);

for i=1:nbrB
    CM(1,1,i) = mean(kol); %Masscentrum x-led 
    CM(2,1,i) = mean(rad); %Masscentrum y-led
end

for k=1:nbrB
    N = size(XY{k},1); % Antal punkter
    XY{k}(:,1) = bsxfun(@minus,XY{k}(:,1),CM(1,1,k)); % Position-CM
    XY{k}(:,2) = bsxfun(@minus,XY{k}(:,2),CM(2,1,k)); % Position-CM
    for i=1:2
        for j=1:2
            G{k}(i,j) = 1/N*dot(XY{k}(:,i),XY{k}(:,j)); %Beräkna gyration tensor
        end
    end
end

%% 

clf
nbrB = size(A,3);%Antal bilder
Trace = zeros(nbrB,1); % Spåret av G
Egenvarde = cell(nbrB,1); % Egenvärde av G
Egenvektorer = cell(nbrB,1); % Egenvektor av G

for i=1:nbrB;
Trace(i,1) = trace(G{i});
[Egenvarde{i},Egenvektorer{i}] = eig(G{i});
end

% Plotta trace som funktion av bilder
figure(1)
plot(Trace) 

% Plotta egenvärde som funktion av bilder
figure(2)
for i=1:nbrB
    plot(i,Egenvarde{i}(1),'*k');hold on;
    plot(i,Egenvarde{i}(2),'vr')
end
legend('Egenvärde 1','Egenvärde 2')



