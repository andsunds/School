%Two-point correlation function

clear all; clf;clc;
filnamn=cell(1,4);
filnamn{1}='data/confined_270304-6-28-min.mat';
filnamn{2}='data/confined_280204-2-32min.mat';
filnamn{3}='data/nonconfined_180304-1-5min.mat';
filnamn{4}='data/nonconfined_250104-1-167min.mat';

A = importdata(filnamn{3});
nbrB = size(A,3);%Antal bilder

XY=cell(nbrB,1); % Cell med positioner

for i=1:nbrB
    [rad,kol] = find(255==A(:,:,i)); % Hitta index f�r position
    XY{i}=[kol,rad];
end

CM = zeros(2,nbrB);

for i=1:nbrB
    CM(1,i) = floor(mean(XY{i}(:,1))); %Masscentrum x-led 
    CM(2,i) = floor(mean(XY{i}(:,2))); %Masscentrum y-led
end
for k=1:nbrB
    XY{k}(:,1) = bsxfun(@minus,XY{k}(:,1),CM(1,k)); % Position relativt CM
    XY{k}(:,2) = bsxfun(@minus,XY{k}(:,2),CM(2,k)); % Position relativt CM
end



s = 50; % Storlek p� unders�kt intervall (z1-z2)
t = 100; % Storlek p� unders�kt intervall (t1-t2)

G = zeros(t+1,2*s+1);
%P = zeros(t+1,2*s+1);

for z = -s:s
    for tau = 0:t
        
        for k = -s:s
            for j = 0:t
            G(k+s+1,j+1) = XY{tau+1}(z+s+1,2)*XY{tau+j+1}(z+2*s+k+1,2)/((2*s+1)*(t+1));
            end
        end        
    end
end
        
%T = linspace(0,t,t+1);
%S = linspace(-s/2,s/2,s+1);
surf(G)

