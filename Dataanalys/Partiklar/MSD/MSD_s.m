function [ s, std_MSD] = MSD_s( fil, N_steps  )
%Beräknar MSD enligt
% s(dt) = 1/(#particles) * sum( ((x(t)-x(0)).^2 ) over all particles

addpath('../')

load('filnamn.mat')
load(['../', kompl], 'intensitet', '-mat')

data =load(['../',filnamn{fil}]);
C = separera(data);






koef=storleksanpassning( fil );

s=zeros(N_steps,1);

index=find(cellfun('length',C)==N_steps).'; 

tmp=zeros(N_steps, length(index));

for j=1:length(index)
    i=index(j);
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    tmp(:,j)=sum(TN.^2,2)/(koef(1)*intensitet{fil}(i).^koef(2));%normera;
    
    %bygger "medelvärde"
    %s=s+tmp/(koef(1)*intensitet{fil}(i).^koef(2));%normerat medelvärde
    
end

%s=s/length(index);
s=sum(tmp,2)/length(index);

std_MSD=std(tmp,0,2)/sqrt(size(tmp,2));


end

