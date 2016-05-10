%% akorr
clc; 
figure(5), clf; 
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



korr=zeros(N,1);%init
korr_sq=zeros(N,1);%init

index=find(cellfun('length',C)==N).'; 
L_i=length(index);
tic

for i=index
    %TN=koordinatbyte(C{i}(1:N,2:3));%positionsdata i TN-koordinater
    TN=(C{i}(1:N,2:3));%positionsdata i TN-koordinater
    
    %Nedan har vi den magiska korrelationsfunktionsberäknaren
    tmpx=triu(TN(:,1)*TN(:,1).');
    tmpy=triu(TN(:,2)*TN(:,2).');
    korr=korr+sum([tmpx(INDEX),tmpy(INDEX)], 2)/(I(i).^koef(2));
    korr_sq=korr_sq+sum([tmpx(INDEX),tmpy(INDEX)].^2, 2)/(I(i).^koef(2));
end
toc


LEN=2*fliplr(1:N).';
korr=korr(1:end)./LEN/L_i;
korr_sq=korr_sq(1:end)./LEN/L_i;

std_korr=sqrt(korr_sq-korr.^2)./sqrt(2*L_i*(LEN-1));

std_korr=std_korr/korr(1);

korr=korr/korr(1);


dt=(0:(N-1)).'*1e-2;

%vilka punkter som ska undersökas
start=2;
stop =700;


c=[ones(stop+1-start, 1) dt(start:stop)]\log(korr(start:stop));%anpassar exponentialsamband

%fprintf('%s\nErhållen tidskonstant: T: %1.3f s,  N: %1.3f s \n', filnamn{fil}, -1/c(1), -1/c(2))
fprintf('Med medelvärdet: %1.3f s \n\n', -1./c(end))


subplot(1,2,fil)

konfidens=0.01;
konfidensfaktor=(norminv(1-konfidens/2,0,1)-norminv(konfidens/2,0,1))/2;

plot(dt, korr, dt,...
     korr+std_korr*konfidensfaktor, '--k', dt, korr-std_korr*konfidensfaktor, '--k')

hold on
plot(dt, exp(c(1)+ c(2)*dt) )
title(filnamn{fil}(6:end-4))
set(gca, 'fontsize', 15, 'yscale', 'log', 'xscale', 'lin')
axis([0,10, .1, 1])
pause(.1)
end


%% simulering
clc; 
figure(5), clf; 
clearvars



N_step=1000;
N_sim=100;

addpath('../');%Lägger till så att create_indecis kan användas
INDEX=create_indecis(N_step);%Tar fram index som sorterar längs med diagonalerna


korr_matrix=zeros(N_step,1);%init


tic
for i=1:N_sim; 
    %TN=koordinatbyte(C{i}(1:N,2:3));%positionsdata i TN-koordinater
    TN=cumsum(randn(N_step,1));
    
    %Nedan har vi den magiska korrelationsfunktionsberäknaren
<<<<<<< HEAD
    tmpx=triu(TN*TN.');
    korr=korr+sum(tmpx(INDEX), 2);
=======
    tmp=triu(TN*TN.');
    korr_matrix=korr_matrix+sum(tmp(INDEX), 2);
>>>>>>> 150342d2733e42e208aa3e854566cc8d78ed7b36
end
toc
s=1;
korr_matrix=korr_matrix(s:end)./fliplr(s:N_step).'/N_sim;


dt=(0:(N_step-1)).'*1e-2;

%vilka punkter som ska undersökas
start=2;
stop =1000;




c=dt(start:stop)\log(korr_matrix(start:stop, :));%anpassar exponentialsamband

%fprintf('%s\nErhållen tidskonstant: T: %1.3f s,  N: %1.3f s \n', filnamn{fil}, -1/c(1), -1/c(2))
fprintf('Med medelvärdet: %1.3f s \n\n', mean(-1./c))



plot(dt(1:stop), korr_matrix(1:stop, :)), hold on
%plot(dt(1:stop), exp( bsxfun(@times,dt(1:stop), c) ) )

set(gca, 'fontsize', 20, 'yscale', 'lin', 'xscale', 'lin')
pause(.1)
