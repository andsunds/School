%% Kovariansundersökning
clc;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom'; 
filnamn{2}='confined_32min_polynom';
filnamn{3}='nonconfined_5min_polynom';
filnamn{4}='nonconfined_167min_polynom';

fil=4;
load(['data/', filnamn{fil}, '.mat'])

framestep=1;%om vi vill undersöka bilder som ligger glesare

N=floor(size(px, 1)/framestep);%antalet bilder som undersöks

S=linspace(0,1,1000);

n = 128;%antalet punkter att kolla korr. i

A     = zeros(n,N);%init
KOV_A = zeros(n);  %init

tic
P=linspace(0,1,n);%punkter längs med strängen

%Beräkna normalvektorer i alla punkter och tider. 
[~, N0]=tangent_normal(PX_mean, PY_mean, P);

Q0=[polyval(PX_mean, P); polyval(PY_mean, P)];%Den undersökt punkten på medelkurvan
Q0=bsxfun(@minus, Q0,...
    mean([polyval(PX_mean, S); polyval(PY_mean, S)], 2) );%tyngdpunkt i origo

for i=1:N %loopa över alla bilder (tid)
    Q1=[polyval(px(i*framestep,:), P); polyval(py(i*framestep,:), P)]; 
    Q1=bsxfun(@minus, Q1,... 
        mean([polyval(px(i*framestep,:), S); polyval(py(i*framestep,:), S)], 2) );%tyngdpunkt i origo
    
    A(:,i)  = diag(N0.'*(Q1-Q0)); %tidsutv. i dim. 2
    KOV_A=KOV_A+A(:,i)*A(:,i).';
end
KOV_A=KOV_A/N;%medelvärde
toc

tic
[V, D] = eig(KOV_A);%tar fram egenvektorer och egenvärden
toc


B=V.'*A;%tidsutveckling längs dim. 2
%B = [B1(1) B1(2), B1(3) ...
%     B2(1) B1(2), B1(3) ...
%     B3(1) B3(2), B3(3) ...
%     ...                   ]
% Där Bi(t) är den i:te moden i bild/tid t.


figure(1), clf
plot(A.'), title('A', 'fontsize', 20)%Massa bröte bara
figure(2), clf
plot(B.'), title('B', 'fontsize', 20)%mindre bröte, men ändå otydligt
%% Fouriertransform av egenvektorerna
clc
figure(3);clf
typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel

index=(n-10):(n-5);
L=arclength(PX_mean, PY_mean);

subplot(2,1,1)%längddomänen
plot(P*L, V(:,index)) 
axis([0, L, min(min(V(:,index)))*1.1, max(max(V(:,index)))*1.1])

typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel
xlabel('$l$ /[px]','Interpreter','Latex');
ylabel('$\mathbf{v}(l)$','Interpreter','Latex')
set(gca,'Fontsize',14, 'yscale','lin', 'xscale', 'lin');


subplot(2,1,2)%frekvensdomänen 

Spektr=fft(V(:,index), [],1);
Fs=n/L;
k=(0:(n/2-1))*Fs/(n-1);

plot(k, (abs(Spektr(1:(n/2), :))));
axis([0, 3e-2, 0, 10])


title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel
xlabel('$k$ /[px$^{-1}$]','Interpreter','Latex');
ylabel('$\mathcal{F}[\mathbf{v}_i](k)$','Interpreter','Latex')
set(gca,'Fontsize',14, 'yscale','lin', 'xscale', 'lin');


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


figure(4), clf
plot(0:(N-1), K) %Häftig korrelatinosfunktioner

typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel

xlabel('$\Delta t$','Interpreter','Latex');
ylabel('$<B_i(t)B_i(t+\Delta t)>_{t}$','Interpreter','Latex')
set(gca,'Fontsize',16)%, 'yscale','log', 'xscale', 'lin');


%% Korskorrelation 
clc;
%Lägger till så att create_indecis kan användas
addpath('../');
%Tar fram index som sorterar längs med diagonalerna
INDEX=create_indecis(N);

i=10;
j=20;

tic
%   /-- bara övre        
%   V    diagonalerna
tmp=triu(B(i,:).'*B(j,:));
%                ^-- korrelationsfunktion på diagonalerna
%                         %v-- medelvärde över diagonalerna
K_kors=sum( tmp(INDEX), 2)./fliplr(1:N).';
fprintf('sum(K_kors)/sqrt(sum(K_kors.^2)) = %.3f\n\n',...
         sum(K_kors)/sqrt(sum(K_kors.^2)))%borde inte denhär vara nära 0?

tmp=triu(B(i,:).'*B(i,:));
K_auto1=sum( tmp(INDEX), 2)./fliplr(1:N).';

tmp=triu(B(j,:).'*B(j,:));
K_auto2=sum( tmp(INDEX), 2)./fliplr(1:N).';
fprintf('sum(K_kors.^2)/max(sum([K_auto1, K_auto2].^2, 1)) = %.3f\n\n',...
         sum(K_kors.^2)/max(sum([K_auto1, K_auto2].^2, 1)))
     
toc

%plottning
figure(5), clf
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

%% Kovarians för B
clc;
KOV_B=zeros(n,n);%init.
tic
for i=1:N %loopa över tid
KOV_B=KOV_B+B(:,i)*B(:,i).';
end
KOV_B=KOV_B/N;%Medelvärde är summan delat på antalaet
toc

%plottning
figure(5); clf
plot(diag(KOV_B))
hold on;grid on
plot(abs(diag(D)))
plot(abs(diag(KOV_B-D)))

xlabel('Sortering', 'interpreter', 'latex')
ylabel('Sorterade egenv\"a{}rden /[px$^2$]', 'interpreter', 'latex')

set(gca, 'fontsize', 20, 'yscale', 'log')

fprintf('\nHur diagonal är KOV_B? \n   Z = %1.2d   (äkta diagonal har Z=0)\n\n', ...
        sqrt(sum(sum((KOV_B-diag(diag(KOV_B))).^2))/sum(diag(KOV_B).^2)) )











