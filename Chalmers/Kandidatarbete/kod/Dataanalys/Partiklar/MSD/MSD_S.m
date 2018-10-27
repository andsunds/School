function [ S, std_MSD ] = MSD_S( fil, N_steps  )
%Beräknar MSD enligt
% S(dt)=(1/T) sum((f(t)-f(t+dt)).^2) over all t

addpath('../')

load('filnamn.mat', 'filnamn', 'kompl')
load(['../', kompl], 'intensitet', '-mat')

data =load(['../',filnamn{fil}]);
C = separera(data);




addpath('../../');%Lägger till så att create_indecis kan användas
INDECIES=create_indecis(N_steps);
LENGHTS=fliplr(1:N_steps).';




koef=storleksanpassning( fil );

S=zeros(N_steps,1);
S_sq=zeros(N_steps,1);

index=find(cellfun('length',C)==N_steps).'; 
L_i=length(index);

%S_matrix=zeros(N_steps, L_i);

for j=1:L_i
    i=index(j);
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    
    tmp=zeros(N_steps,2);
    
    %Hög minnesåtgång
    A=repmat(TN(:,1), 1,N_steps);
    A=triu(A.'-A);
    tmp(:,1)=sum(A(INDECIES).^2, 2);
    
    A=repmat(TN(:,2), 1,N_steps);
    A=triu(A.'-A);
    tmp(:,2)=sum(A(INDECIES).^2, 2);
    
    tmp=tmp/(koef(1)*intensitet{fil}(i).^koef(2));
    
    S=S+sum(tmp,2);
    S_sq=S_sq+sum(tmp.^2,2);
    
    %S_matrix(:,j)=sum(tmp,2)/(koef(1)*intensitet{fil}(i).^koef(2));
end

%S=sum(S_matrix,2)/L_i;
%std_S=std(S_matrix,0, 2)/sqrt(L_i);

S=S./(L_i*LENGHTS);
S_sq=S_sq./(L_i*LENGHTS);

std_MSD=sqrt(S_sq-S.^2)./sqrt(L_i*(LENGHTS-1));

end
