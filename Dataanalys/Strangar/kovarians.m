%% Kovariansundersökning
clc;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom'; 
filnamn{2}='confined_32min_polynom';
filnamn{3}='nonconfined_5min_polynom';
filnamn{4}='nonconfined_167min_polynom';

fil=4;
load(['data/', filnamn{fil}, '.mat'])

N=size(px, 1);
S=linspace(0,1,1000);

n=100;%antalet punkter att kolla korr. i

A=zeros(n,N);%init
M=zeros(n);

tic
P=linspace(0,1,n);%loopa över hela strängen

%Beräkna normalvektorer i alla punkter och tider. 
[~, N0]=tangent_normal(PX_mean, PY_mean, P);

Q0=[polyval(PX_mean, P); polyval(PY_mean, P)];%Den undersökt punkten på medelkurvan
Q0=bsxfun(@minus, Q0,...
    mean([polyval(PX_mean, S); polyval(PY_mean, S)], 2) );%tyngdpunkt i origo

for i=1:N %loopa över alla bilder (tid)
    Q1=[polyval(px(i,:), P); polyval(py(i,:), P)]; 
    Q1=bsxfun(@minus, Q1,... 
        mean([polyval(px(i,:), S); polyval(py(i,:), S)], 2) );%tyngdpunkt i origo
    
    A(:,i)  = diag(N0.'*(Q1-Q0)); %tidsutv. i dim. 2
    M=M+A(:,i)*A(:,i).';
end
M=M/N;%medelvärde
toc

tic
[V, ~] = eig(M);%tar fram egenvektorer
toc


B=V.'*A;
%B = [B1(1) B1(2), B1(3) ...
%     B2(1) B1(2), B1(3) ...
%     B3(1) B3(2), B3(3) ...
%     ...                   ]
% Där Bi(t) är den i:te moden i bild/tid t.


%figure(1), clf
%plot(A.') %Massa bröte bara
%figure(2), clf
%plot(B.') %mindre bröte, men ändå otydligt

%% Autokorrelation
clc
K=zeros(N,n);%init.
%Lägger till så att create_indecis kan användas
addpath('../');
%Tar fram index som sorterar längs med diagonalerna
INDEX=create_indecis(N);

tic
for i=1:n;
%   /-- bara övre        
%   V    diagonalerna
tmp=triu(B(i,:).'*B(i,:));
%                ^-- korrelationsfunktion på diagonalerna
%                         %v-- medelvärde över diagonalerna
K(:,i)=sum( tmp(INDEX), 2)./fliplr(1:N).';
end
toc


figure(3), clf
plot(0:(N-1), K) %Häftig korrelatinosfunktioner

typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel

xlabel('$\Delta t$','Interpreter','Latex');
ylabel('$<B_i(t)B_i(t+\Delta t)>_{t}$','Interpreter','Latex')
set(gca,'Fontsize',16)%, 'yscale','log');


%% Korskorrelation
clc;
K=zeros(N,1);%init.
%Lägger till så att create_indecis kan användas
addpath('../');
%Tar fram index som sorterar längs med diagonalerna
INDEX=create_indecis(N);

i=60;
j=50;

tic
%   /-- bara övre        
%   V    diagonalerna
tmp=triu(B(i,:).'*B(j,:));
%                ^-- korrelationsfunktion på diagonalerna
%                         %v-- medelvärde över diagonalerna
K_kors=sum( tmp(INDEX), 2)./fliplr(1:N).';
fprintf('sum(K_kors)/sqrt(sum(K_kors.^2)) = %.3f\n\n',...
         sum(K_kors)/sqrt(sum(K_kors.^2)))%borde inte denhär vara nära 0?
fprintf('sum(K_kors.^2)/max(sum([K_auto1, K_auto2].^2, 1)) = %.3f\n\n',...
         sum(K_kors.^2)/max(sum([K_auto1, K_auto2].^2, 1)))
     
tmp=triu(B(i,:).'*B(i,:));
K_auto1=sum( tmp(INDEX), 2)./fliplr(1:N).';

tmp=triu(B(j,:).'*B(j,:));
K_auto2=sum( tmp(INDEX), 2)./fliplr(1:N).';
toc


figure(4), clf
plot(0:(N-1), K_kors) 
hold on
plot(0:(N-1), K_auto1) 
plot(0:(N-1), K_auto2) 

typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel

leg=legend(sprintf('$i=%d$, $j=%d$', i, j), ...
           sprintf('$i=%d$, $j=%d$', i, i), ...
           sprintf('$i=%d$, $j=%d$', j, j));
set(leg, 'interpreter', 'Latex')

xlabel('$\Delta t$','Interpreter','Latex');
ylabel('$<B_i(t)B_j(t+\Delta t)>_{t}$','Interpreter','Latex')
set(gca,'Fontsize',16)%, 'yscale','log');

