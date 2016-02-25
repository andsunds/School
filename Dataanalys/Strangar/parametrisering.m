%% Parametrisering av strängar som funktion av båglängd s

clear all; clf;clc;
filnamn=cell(1,4);
filnamn{1}='confined_270304-6-28-min.mat';
filnamn{2}='confined_280204-2-32min.mat';
filnamn{3}='nonconfined_180304-1-5min.mat';
filnamn{4}='nonconfined_250104-1-167min.mat';

A = importdata(filnamn{3});
nbrB = size(A,3);%Antal bilder

XY=cell(nbrB,1); % Cell med positioner
s = cell(nbrB,1); % Båglängd, används till parametrisering s.a 
                  % x=x(s) och y=y(s)
S = cell(nbrB,1); % Normerad variant av s. 
X = cell(nbrB,1); % Spline x(s)
Y = cell(nbrB,1); % Spline y(s)
L = zeros(nbrB,1); % Total längd L
for i=1:nbrB
    [rad,kol] = find(255==A(:,:,i)); % Hitta index för position
    XY{i}=[kol,rad];
    
    % Skapa s-parameter som kumulativ summa av vägskillnader
    s{i} = cumsum(sqrt([0,diff(kol')].^2 + [0,diff(rad')].^2)); 
    
    L(i) = s{i}(end); 
    
    % Normerad variant för att möjliggöra jämförelse % mellan tider av 
    % olika antal punkter och olika totala längder L.
    S{i} = s{i}/max(s{i}); 
                           
                           
    % Ta fram spline x(s) och y(s). Byt ut s{i}->S{i} för normerad variant, dvs
    % att man kan jämföra strängar av olika L. 
    X{i} = spline(S{i},kol); 
    Y{i} = spline(S{i},rad);
end 


%% 
Xs = cell(nbrB,1); % x(s) från spline
Ys = cell(nbrB,1); % y(s) från spline
for i = 1:nbrB
    % ppval tar en spline X{i} och en parametervektor s{i} och beräknar
    % motsvarande x(s)
    Xt{i} = ppval(X{i},s{i}); % x(s) från spline
    Yt{i} = ppval(Y{i},s{i}); % y(s) från spline
    plot(Xt{i},Yt{i})
    pause(0.5)
end



%% Surf
%Välj S{i} vid spline anpassning
clf;
Ss = linspace(0,1,nbrB); % s-parameter 
T = linspace(0,nbrB-1,nbrB); % Tidsvektor   Ss och T behöver va lika stora.
XYT = zeros(nbrB,nbrB,2); % XY-värden som "funktion" av tiden T

for i=1:nbrB
    XYT(i,:,1) = ppval(X{i},Ss); % Beräkna x(Ss)
    XYT(i,:,2) = ppval(Y{i},Ss); % Beräkna y(Ss)
    XYT(i,:,3) = T; % Tidsvektor 
end

surf(XYT(:,:,1),XYT(:,:,2),XYT(:,:,3)) % '3D' variant
xlabel('x')
ylabel('y')
zlabel('t')
set(gca,'Fontsize',30);





