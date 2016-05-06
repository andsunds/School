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
s_sq=zeros(N_steps,1);

index=find(cellfun('length',C)==N_steps).'; 
L_i=length(index);

%s_matrix=zeros(N_steps, L_i);

for j=1:L_i
    i=index(j);
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    
    %s_matrix(:,j)=sum(TN.^2,2)/(koef(1)*intensitet{fil}(i).^koef(2));%normera;
    
    tmp=sum(TN.^2,2)/(koef(1)*intensitet{fil}(i).^koef(2));
    %bygger "medelvärde"
    s=s+tmp;
    s_sq=s_sq+tmp.^2;
end

%s=sum(s_matrix,2)/L_i;
%std_MSD=std(s_matrix,0,2)/sqrt(L_i);

s=s/L_i;

std_MSD=sqrt(s_sq/L_i - s.^2)/sqrt(L_i);


end

