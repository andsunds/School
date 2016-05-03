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



t=(0:(N-1)).'/framerate;

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
plot(t,A.'*1e6), title('A', 'fontsize', 20)%Massa bröte bara
xlabel('$t$ /[s]','Interpreter','Latex');
ylabel('$\mathbf{v}(l) /[\mu \mathrm{m}]$','Interpreter','Latex')
set(gca,'Fontsize',14, 'yscale','lin', 'xscale', 'lin');

figure(2), clf
plot(t,B.'*1e6), title('B', 'fontsize', 20)%mindre bröte, men ändå otydligt
xlabel('$t$ /[s]','Interpreter','Latex');
ylabel('$\mathbf{v}(l) /[\mu \mathrm{m}]$','Interpreter','Latex')
set(gca,'Fontsize',14, 'yscale','lin', 'xscale', 'lin');
%% Fouriertransform av egenvektorerna
clc
figure(3);clf
typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel

index=(n-3):1:(n);

L=arclength(PX_mean, PY_mean);

subplot(2,2,[1,2])%längddomänen
plot(P*L, V(:,index)) 
axis([0, L, min(min(V(:,index)))*1.1, max(max(V(:,index)))*1.1])

typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel
xlabel('$l /[\mathrm{m}]$','Interpreter','Latex');
ylabel('$\mathbf{v}(l)$','Interpreter','Latex')
set(gca,'Fontsize',14, 'yscale','lin', 'xscale', 'lin');
grid on

subplot(2,2,3)%frekvensdomänen 

Spektr=fft(V, [],1);
Fs=n/L;
k0=(0:(n/2-1))*Fs/(n-1);

plot(repmat(k0.', 1, length(index)), (abs(Spektr(1:(n/2), index))));
axis([0, 1e6, 0, 10])

title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel
xlabel('$k /[\mathrm{m}^{-1}]$','Interpreter','Latex');
ylabel('$\mathcal{F}[\mathbf{v}_i](k)$','Interpreter','Latex')
set(gca,'Fontsize',14, 'yscale','lin', 'xscale', 'lin');


subplot(2,2,4)%frekvensberoende

[~, i]=max(abs(Spektr(2:(n/2), :)));

f_max=k0(i);

d=diag(D);

plot(f_max(index), d(index), '*')

title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel
xlabel('$k /[\mathrm{m}^{-1}]$','Interpreter','Latex');
ylabel('var[$B_i] /[\mathrm{m}^2$]','Interpreter','Latex')
set(gca,'Fontsize',14, 'yscale','log', 'xscale', 'lin');

%% Anpassning av vågtal
clc;
figure(3);clf
typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel

index=(n-5):1:(n-0);
%index=n-9;

L=arclength(PX_mean, PY_mean);

subplot(2,1,1)%längddomänen
plot(P*L, V(:,index) ) 
%axis([0, L, min(min(V(:,index)))*1.1, max(max(V(:,index)))*1.1])

typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel
xlabel('$l /[\mathrm{m}]$','Interpreter','Latex');
ylabel('$\mathbf{v}(l)$','Interpreter','Latex')
set(gca,'Fontsize',14, 'yscale','lin', 'xscale', 'lin');
grid on



hold on
Fs=n/L;
k0=(0:(n/2-1))*Fs/(n-1);

k=zeros(length(index),1);

for i=1:length(index)
    j=index(i);
    f=@(X) X(1)+ X(2)*cos( X(3)*P*L+X(4) );
    
    Spektr=fft(V(:,j)-mean(V(:,j)), [],1);
    [~, inx]=max(abs(Spektr(2:(n/2), :)));
    
    if 10*k0(inx)>2e5
        start=[0, 1, 10*k0(inx) ,0];
    else
        start=[0, 1, 1e5 ,0];
    end
    [K, fval]=fminsearch(@(X) sum((f(X)-V(:,j).').^2) , start, ...
        optimset('MaxFunEvals',4000, 'MaxIter',1000, 'tolfun', 1e-9));
    
    tol=0.5;
    if fval>tol
        fprintf('Miss på första försöket, \t \t n-index = %d\n',n-j)
    counter=0;
    for korr=linspace(-.5,5,50);
       start(3)=start(3)*(1+korr);
       [K, fval]=fminsearch(@(X) sum((f(X)-V(:,j).').^2) , start, ...
       optimset('MaxFunEvals',4000, 'MaxIter',1000, 'tolfun', 1e-9));
       counter=counter+1;    
       if fval<tol
           break
       end
    end
    fprintf('Använde, %d försök\n',counter)
    end
    
    fprintf('fval = %f,\t n-index = %d\n\n',fval,n-j)
    k(i)=K(3);
    
    plot(P*L, f(K),'--')
end

subplot(2,1,2)%frekvensdomänen 
d=diag(D);

plot(k, d(index), '*')

title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel
xlabel('$k /[\mathrm{m}^{-1}]$','Interpreter','Latex');
ylabel('var[$B_i] /[\mathrm{m}^2$]','Interpreter','Latex')
set(gca,'Fontsize',14, 'yscale','log', 'xscale', 'lin');


%% Autokorrelation
clc;clf;
K=zeros(N,n);%init.
%Lägger till så att create_indecis kan användas
addpath('../');
%Tar fram index som sorterar längs med diagonalerna
INDEX=create_indecis(N);

dt=t;

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
plot(dt, K) %Häftig korrelatinosfunktioner

typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel

xlabel('$\Delta t /[\mathrm{s}]$','Interpreter','Latex');
ylabel('$<B_i(t)B_i(t+\Delta t)>_{t} /[\mathrm{m}^2]$','Interpreter','Latex')
set(gca,'Fontsize',16, 'yscale','log', 'xscale', 'lin');

%%  Anpassning av relaxationstider
% K�r alla delar ovan innan denna.
% h best�mmer antalet tidssteg f�r mod 1,2,3,4
t=6;
h1 = [9 6 22 4 3 0 2 2 3].*[repmat(1,1,t) repmat(0,1,9-t)];%Fil1
h2 = [4 12 4 3];%Fil2
h3 = [15 21 14 6 4 4 3 2 10];%.*[repmat(1,1,t) repmat(0,1,9-t)];%Fil3
h4 = [35 15 12 4 2]; %Fil4 
if fil==1
    h=h1;
    elseif fil==2
    h=h2;
    elseif fil==3
    h=h3;
    else
    h=h4;
end


k0=size(h,2); % Antal moder som vill unders�kas
Y = zeros(N,k0);
Index = (n-k0+1):(n); % Vilka moder studeras

for j=1:k0
    Y(:,j) = [ones(h(j),1);NaN(N-h(j),1)]; % Placera ettor f�r de intervall som ska va kvar och NaN p� de som ej ska va med.
end
T=fliplr(K(:,Index)).*Y; % Flippa K f�r att placera f�rsta moden i kolumn 1, andra i kolumn 2 etc.
                         % Spara tillr�ckligt m�nga element dt

figure(5), clf
subplot(2,2,[1,2])
plot(dt, T) %Häftig korrelatinosfunktioner

xlabel('$\Delta t /[\mathrm{s}]$','Interpreter','Latex');
ylabel('$<B_i(t)B_i(t+\Delta t)>_{t} /[\mathrm{m}^2]$','Interpreter','Latex')
set(gca,'Fontsize',16, 'yscale','log', 'xscale', 'lin');
hold on;
tau = zeros(1,k0);
nbr = zeros(1,k0);

for i=1:k0
    x = dt(1:h(i));
    y = log(T(1:h(i),i));
    p = polyfit(x,y,1);
    
    tau(i) = abs(1/p(1)); % Ordnad s.a. tau(1) motsvarar mod med st�rst egenv�rde
    D = exp(p(1).*x).*exp(p(2));
    plot(dt(1:h(i)),D,'--')
end

subplot(2,2,3)
vagtal = fliplr(f_max(Index));
plot(vagtal,tau,'*r')
axis([0 10^5 0 15])

set(gca,'yscale','log','xscale','lin','Fontsize',14)
ylabel('$\tau_n$','Interpreter','Latex','Fontsize',26)
xlabel('$k /[\mathrm{m}^{-1}]$','Interpreter','Latex');

subplot(2,2,4)
plot(f_max(Index), d(Index), '*')
xlabel('$k /[\mathrm{m}^{-1}]$','Interpreter','Latex');
ylabel('var[$B_i] /[\mathrm{m}^2$]','Interpreter','Latex')
set(gca,'Fontsize',14, 'yscale','log', 'xscale', 'lin');



disp(fliplr(f_max(index)))
%save('relaxtiderfil4.mat','vagtal','tau')


figure(6)
plot(fliplr(d(Index)'),tau,'*')
set(gca,'yscale','lin','Fontsize',26)
ylabel('$\tau_n /[\mathrm{s}]$','Interpreter','Latex');
xlabel('var[$B_i] /[\mathrm{m}^2$]','Interpreter','Latex')
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
figure(7), clf
plot(dt, K_kors) 
hold on
plot(dt, K_auto1) 
plot(dt, K_auto2) 

typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel

leg=legend(sprintf('$i=%d$, $j=%d$', i, j), ...
           sprintf('$i=%d$, $j=%d$', i, i), ...
           sprintf('$i=%d$, $j=%d$', j, j));
set(leg, 'interpreter', 'Latex')

xlabel('$\Delta t /[\mathrm{s}]$','Interpreter','Latex');
ylabel('$<B_i(t)B_i(t+\Delta t)>_{t} /[\mathrm{m}^2]$','Interpreter','Latex')
set(gca,'Fontsize',16)%, 'yscale','log', 'xscale', 'lin');


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
figure(8); clf
plot(diag(KOV_B))
hold on;grid on
plot(abs(diag(D)))
plot(abs(diag(KOV_B-D)))

xlabel('Sortering', 'interpreter', 'latex')
ylabel('Sorterade egenv\"a{}rden /[$\mathrm{m}^2$]', 'interpreter', 'latex')

set(gca, 'fontsize', 20, 'yscale', 'log')

fprintf('\nHur diagonal är KOV_B? \n   Z = %1.2d   (äkta diagonal har Z=0)\n\n', ...
        sqrt(sum(sum((KOV_B-diag(diag(KOV_B))).^2))/sum(diag(KOV_B).^2)) )











