%% akorr
clc; 
figure(3), clf; 
clearvars
load('filnamn.mat')
load(kompl, '-mat')

fil=1;
N=1000;


addpath('../');%Lägger till så att create_indecis kan användas
INDEX=create_indecis(N);%Tar fram index som sorterar längs med diagonalerna

for fil=1:2 %för att undersöka datan för båda celltyperna
%laddar in data
data =load(filnamn{fil});
C = separera(data);

I=intensitet{fil};

%Beräklna fram en normering för variansen som fkn av intensitet
lambda=sqrt(std_n{fil}.^2+std_t{fil}.^2 -2*sigma_brus(fil).^2 );
%beräkna koefficienter för lutning i loglog
koef=[ones(size(I)), log(I)]\log(lambda);
koef(1)=exp(koef(1));%konverterar till potenssamband


korr_T=zeros(N,1);%init
korr_N=zeros(N,1);%init

tic
for i=find(cellfun('length',C)==N).'; 
    %TN=koordinatbyte(C{i}(1:N,2:3));%positionsdata i TN-koordinater
    TN=(C{i}(1:N,2:3));%positionsdata i TN-koordinater
    
    %Nedan har vi den magiska korrelationsfunktionsberäknaren
    tmp1=triu(TN(:,1)*TN(:,1).');
    tmp2=sum(tmp1(INDEX), 2);
    korr_T=korr_T+tmp2/(I(i).^koef(2));
    
    tmp1=triu(TN(:,2)*TN(:,2).');
    tmp2=sum(tmp1(INDEX), 2);
    korr_N=korr_N+tmp2/(I(i).^koef(2));
end
toc
korr_T=korr_T./fliplr(1:N).';
korr_N=korr_N./fliplr(1:N).';

korr=[korr_T/korr_T(1), korr_N/korr_N(1)];

dt=(0:(N-1)).'*1e-2;

%vilka punkter som ska undersökas
start=2;
stop =200;

c=dt(start:stop)\log(korr(start:stop, :));%anpassar exponentialsamband

fprintf('%s\nErhållen tidskonstant: T: %1.3f s,  N: %1.3f s \n', filnamn{fil}, -1/c(1), -1/c(2))
fprintf('Med medelvärdet: %1.3f s \n\n', mean(-1./c))


subplot(1,2,fil)
plot(dt(1:stop), korr(1:stop, :)), hold on
plot(dt(1:stop), exp( bsxfun(@times,dt(1:stop), c))  )
set(gca, 'fontsize', 20, 'yscale', 'log', 'xscale', 'lin')
pause(.1)
end
