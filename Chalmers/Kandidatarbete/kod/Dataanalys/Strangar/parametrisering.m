%% Parametrisering av str�ngar som funktion av b�gl�ngd s

clear all; clf;clc;
filnamn=cell(1,4);
filnamn{1}='data/confined_270304-6-28-min.mat';
filnamn{2}='data/confined_280204-2-32min.mat';
filnamn{3}='data/nonconfined_180304-1-5min.mat';
filnamn{4}='data/nonconfined_250104-1-167min.mat';

A = importdata(filnamn{3});
nbrB = size(A,3);%Antal bilder

XY=cell(nbrB,1); % Cell med positioner
s = cell(nbrB,1); % B�gl�ngd, anv�nds till parametrisering s.a 
                  % x=x(s) och y=y(s)
S = cell(nbrB,1); % Normerad variant av s. 
X = cell(nbrB,1); % Spline x(s)
Y = cell(nbrB,1); % Spline y(s)
L = zeros(nbrB,1); % Total l�ngd L
for i=1:nbrB
    [rad,kol] = find(255==A(:,:,i)); % Hitta index f�r position
    XY{i}=[kol,rad];
    
    % Skapa s-parameter som kumulativ summa av v�gskillnader
    s{i} = cumsum(sqrt([0,diff(kol')].^2 + [0,diff(rad')].^2)); 
    
    L(i) = s{i}(end); 
    
    % Normerad variant f�r att m�jligg�ra j�mf�relse % mellan tider av 
    % olika antal punkter och olika totala l�ngder L.
    S{i} = s{i}/max(s{i}); 
                           
                           
    % Ta fram spline x(s) och y(s). Byt ut s{i}->S{i} f�r normerad variant, dvs
    % att man kan j�mf�ra str�ngar av olika L. 
    X{i} = spline(S{i},kol); 
    Y{i} = spline(S{i},rad);
end 

clear A;

%% 
Xs = cell(nbrB,1); % x(s) fr�n spline
Ys = cell(nbrB,1); % y(s) fr�n spline
for i = 1:nbrB
    % ppval tar en spline X{i} och en parametervektor s{i} och ber�knar
    % motsvarande x(s)
    SS=linspace(0,1,10000);
    Xt{i} = ppval(X{i},SS); % x(s) fr�n spline
    Yt{i} = ppval(Y{i},SS); % y(s) fr�n spline
    plot(Xt{i},Yt{i})
    pause(0.5)
end



%% Surf
%V�lj S{i} vid spline anpassning
clf;
Ss = linspace(0,1,nbrB); % s-parameter 
T = linspace(0,nbrB-1,nbrB); % Tidsvektor   Ss och T beh�ver va lika stora.
XYT = zeros(nbrB,nbrB,2); % XY-v�rden som "funktion" av tiden T

for i=1:nbrB
    XYT(i,:,1) = ppval(X{i},Ss); % Ber�kna x(Ss)
    XYT(i,:,2) = ppval(Y{i},Ss); % Ber�kna y(Ss)
    XYT(i,:,3) = T; % Tidsvektor 
end

surf(XYT(:,:,1),XYT(:,:,2),XYT(:,:,3)) % '3D' variant
xlabel('x')
ylabel('y')
zlabel('t')
set(gca,'Fontsize',30);





