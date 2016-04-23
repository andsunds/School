function [ S ] = MSD_S( fil, N_steps  )
%Beräknar MSD enligt
% S(dt)=(1/T) sum((f(t)-f(t+dt)).^2) over all t


load('../kompleterande_data.mat', 'filnamn', 'intensitet', '-mat')


addpath('../../');%Lägger till så att create_indecis kan användas
INDECIES=create_indecis(N_steps);
LENGHTS=fliplr(1:N_steps).';


data =load(['../',filnamn{fil}]);
C = separera(data);



koef=storleksanpassning( fil );

S=zeros(N_steps,1);

index=find(cellfun('length',C)==N_steps).'; 
for i=index
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    
    tmp=zeros(N_steps,2);
    
    %Hög minnesåtgång
    A=repmat(TN(:,1), 1,N_steps);
    A=triu(A.'-A);
    tmp(:,1)=sum(A(INDECIES).^2, 2)./LENGHTS;
    
    A=repmat(TN(:,2), 1,N_steps);
    A=triu(A.'-A);
    tmp(:,2)=sum(A(INDECIES).^2, 2)./LENGHTS;
    
    
    % Detta är samma normering som för rörligheten, 
    % kanske skulle man hitta på något annat.
    S=S+sum(tmp,2)/(koef(1)*intensitet{fil}(i).^koef(2));
end

S=S/length(index);

end

