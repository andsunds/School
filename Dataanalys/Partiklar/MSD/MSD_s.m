function [ s ] = MSD_s( fil, N_steps )
%Beräknar MSD enligt
% s(dt) = 1/(#particles) * sum( ((x(t)-x(0)).^2 ) over all particles

addpath('../')

load('filnamn.mat')
load(['../', kompl], 'intensitet', '-mat')

data =load(['../',filnamn{fil}]);
C = separera(data);


koef=storleksanpassning( fil );

%initialisering
s=zeros(N_steps,1);

%loop över alla partiklar med tidsdata upp till 10s
index=find(cellfun('length',C)==N_steps).'; 
for i=index;
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    
    %bygger "medelvärde"
    s=s+sum(TN.^2, 2)/(koef(1)*intensitet{fil}(i).^koef(2));%normerat medelvärde
    
end

s=s/length(index);



end

