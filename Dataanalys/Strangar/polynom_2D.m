%% 2D-polynom
clc;clf;clear all

%Ladda in data
filnamn=cell(1,4);
filnamn{1}='confined_28min.mat'; 
filnamn{2}='confined_32min.mat';
filnamn{3}='nonconfined_5min.mat';
filnamn{4}='nonconfined_167min.mat';

load(['data/' filnamn{4}]);



N=size(coordinates, 1);%antalet bilder

%gradtal för polynomanpassning
grad=20;
px=zeros(N,grad+1);%init
py=zeros(N,grad+1);%init

ws=warning('off', 'all');%Stänger av polyfits varingar
tic
for i=1:N
%Kurvlängds vekot med varje punkts pos. längs kurvan.
s=L_string(i,1:N_points(i))/L_string(i,N_points(i));
%%%%%%%%                   ^alla blir lika långa

x=coordinates(i,1:N_points(i), 1);
y=coordinates(i,1:N_points(i), 2);
%polynomanpassar x och y efter kurvlängden
px(i,:) = polyfit(s,x,grad);
py(i,:) = polyfit(s,y,grad);

% %plottning:
% S=linspace(0, 1,10000);
% X=polyval(px(i,:),S);
% Y=polyval(py(i,:),S);
% 
plot(x-mean(x),y-mean(y)), hold on
% plot(X,Y)
% hold off
% pause(.5)
end
toc
warning(ws)

% Medelvärde av polynom är samma som medevärde av koef.
PX=mean(px,1);
PY=mean(py,1);
%Tarfram x- och y-värden som fkn av kurvlängd
S=linspace(0, 1 ,1000);
X=polyval(PX,S);
Y=polyval(PY,S);
%plotta
plot(X-mean(X),Y-mean(Y), 'K', 'linewidth',10)
axis equal

%% Film med sträng och medelpositionen
% Alla bilder är centrerade kring medelvärdet
clc;clf

plot(X-mean(X),Y-mean(Y), 'k', 'linewidth',4)%plottar medelvärdet
hold on
%Första steget i for-loopen är lite speciellt
i=1;
XP=polyval(px(i,:),S);
YP=polyval(py(i,:),S);
x=coordinates(i,1:N_points(i), 1);
y=coordinates(i,1:N_points(i), 2);

h=plot(x-mean(XP),y-mean(YP));
p=plot(XP-mean(XP),YP-mean(YP));
pause(.1)
%sen fortsätter filmen
for i=2:N
    %Uppdatera x- och y-koordinaterna:
    x=coordinates(i,1:N_points(i), 1);
    y=coordinates(i,1:N_points(i), 2);
    %uppdatera polynomanpassningen:
    XP(i,:)=polyval(px(i,:),S);
    YP(i,:)=polyval(py(i,:),S);
    %Snabbare med set(...) när man vill uppdatera bilder.
    set(p, 'XData',XP(i,:)-mean(x) ,'YData',YP(i,:)-mean(y));
    set(h, 'XData',x-mean(x),'YData',y-mean(y));
    %formatering av bilden:
    axis equal, axis([-300,300, -200, 200])
    %pause(.1)
end
hold off

%% avstånd från jmvkt.
n=50;
K=zeros(1,n-1);
Q=500;
avs(:,1)=XP(:,Q)-X(Q);
avs(:,2)=YP(:,Q)-Y(Q);
avs=avs';

for i=1:N;
%Tangenten ges av derivatan:
dX=polyder(X);%derivera x
dY=polyder(Y);%derivera y

jmN=[dX(Q),dY(Q)]/sqrt([dX(Q),dY(Q)].^2);

%Tangentvektor i de specifika pkt.
%T=[polyval(dx, l); polyval(dy, l) ];
avs=avs./repmat(sqrt(sum(avs.^2,1)),2,1);%normering

%Beräkna korr
for s=1:n
for dl=1:(n-s)
    %Sumerar delar till medelvärden av korrelationsfunktionen
    K(dl)=K(dl)+(avs(:,s).'*avs(:,s+dl-1))/(n-dl)/N;
end
end

end

%% Tangentvektorskorr (Finns förbättringspotential)
%<t(s) * t(s+l)> ~ exp(-l/L_P)
clc;clf

n=100;%antalet punkter att kolla korr. i

K=zeros(1,n-1);%init.
l=linspace(0,1,n);%Vilka punkter vi ska kolla efter tangentvektor

tic
for i=1:N;
%Tangenten ges av derivatan:
dx=polyder(px(i,:));%derivera x
dy=polyder(py(i,:));%derivera y

%Tangentvektor i de specifika pkt.
T=[polyval(dx, l); polyval(dy, l) ];
T=T./repmat(sqrt(sum(T.^2,1)),2,1);%normering

%Beräkna korr
for s=1:n
for dl=1:(n-s)
    %Sumerar delar till medelvärden av korrelationsfunktionen
    K(dl)=K(dl)+(T(:,s).'*T(:,s+dl-1))/(n-dl)/N;
end
end

end
toc
%plottar
plot(linspace(0,1,n-1), K)

%axis([0,1, 0,1])

%set(gca, 'fontsize',15, 'yscale', 'log')































