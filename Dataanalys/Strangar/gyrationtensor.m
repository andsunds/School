% Gyration tensor
clear all; clf;clc;
filnamn=cell(1,4);
filnamn{1}='data/confined_270304-6-28-min.mat';
filnamn{2}='data/confined_280204-2-32min.mat';
filnamn{3}='data/nonconfined_180304-1-5min.mat';
filnamn{4}='data/nonconfined_250104-1-167min.mat';

A = importdata(filnamn{3});
nbrB = size(A,3);%Antal bilder

XY=cell(nbrB,1); % Cell med positioner
G = cell(nbrB,1); % Gyration tensor

for i=1:nbrB
    [rad,kol] = find(255==A(:,:,i)); % Hitta index f�r position
    XY{i}=[kol,rad];
end

CM = zeros(2,nbrB);

for i=1:nbrB
    CM(1,i) = mean(XY{i}(:,1)); %Masscentrum x-led 
    CM(2,i) = mean(XY{i}(:,2)); %Masscentrum y-led
end

for k=1:nbrB
    N = size(XY{k},1); % Antal punkter
    XY{k}(:,1) = bsxfun(@minus,XY{k}(:,1),CM(1,k)); % Position relativt CM
    XY{k}(:,2) = bsxfun(@minus,XY{k}(:,2),CM(2,k)); % Position relativt CM
    for i=1:2
        for j=1:2
            G{k}(i,j) = 1/N*dot(XY{k}(:,i),XY{k}(:,j)); %Ber�kna gyration tensor
        end
    end
end

% Ber�kning av sp�ret av G samt dess egenv�rden
%1:a delen beh�ver k�ras f�rst

clf
nbrB = size(A,3);%Antal bilder
Trace = zeros(nbrB,1); % Sp�ret av G
Egenvarde = cell(nbrB,1); % Egenv�rde av G
Egenvektorer = cell(nbrB,1); % Egenvektor av G

for i=1:nbrB;
Trace(i,1) = trace(G{i});
[Egenvektorer{i},Egenvarde{i}] = eig(G{i});
end

% Plotta trace som funktion av bilder(tid)
figure(1)
plot(Trace) 
title('Tr(G)','Interpreter','Latex')
xlabel('Tid','Interpreter','Latex')
set(gca,'Fontsize',30)

% Plotta egenv�rde som funktion av bilder(tid)
figure(2)
for i=1:nbrB
    plot(i,Egenvarde{i}(1,1),'*k');hold on;
    plot(i,Egenvarde{i}(2,2),'vr')
end
leg = legend('$\lambda_x$','$\lambda_y$');
set(leg,'Interpreter','Latex')
title('Egenv�rden till G','Interpreter','Latex')
xlabel('Tid','Interpreter','latex')
set(gca,'Fontsize',30)


%% Asphericity
%1:a & 2:a delen beh�ver k�ras f�rst

d = 2; % Dimension
lambdaMedel = bsxfun(@times,Trace,1/d); % Medel egenv�rde

I = eye(2);
Ad = zeros(nbrB,1); % Asphericity; =0 f�r sf�riskt symm, =1 f�r rak linje
Ghat = cell(nbrB,1); 
for i=1:nbrB
    Ghat{i} = G{i}-lambdaMedel(i,1)*I;
    Ad(i,1) = d/(d-1)*trace(Ghat{i}^2)/(Trace(i)^2);
end

figure(3)
plot(Ad)
title('Asphericity A$_d = \frac{d}{d-1}\frac{Tr\hat{G}^2}{(TrG)^2}$','Interpreter','Latex')
xlabel('Tid','Interpreter','latex')
set(gca,'Fontsize',30)



% "film" p� str�ngen f�r att kunna j�mf�ra A_d och r�relsen
bilder_s=3; %bilder/sekund
disp('speltid=')
disp(size(A,3)/bilder_s)

%%
figure(4)
for i = 1:size(A,3)
  image(A(:,:,i))
  pause(bilder_s^-1)
  i % Full�sning f�r att se vart i tiden filmen �r
end

%% Fouriertransform av egenv�rden
L = zeros(2,nbrB);

for i=1:nbrB 
    L(1,i) = Egenvarde{i}(1,1); % Lambda_x
    L(2,i) = Egenvarde{i}(2,2); % Lambda_y
end

LX = fft(L(1,:));
LY = fft(L(2,:));


figure(5)
plot(abs(LX).^2)
figure(6)
plot(abs(LY).^2)


