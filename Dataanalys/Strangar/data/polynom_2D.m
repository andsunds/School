%% 2D-polynom
clc;clf;clear all

%Ladda in data
filnamn=cell(1,4);
filnamn{1}='confined_28min.mat'; 
filnamn{2}='confined_32min.mat';
filnamn{3}='nonconfined_5min.mat';
filnamn{4}='nonconfined_167min.mat';

fil=1;
load([filnamn{fil}]);



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
PX_mean=mean(px,1);
PY_mean=mean(py,1);
%Tarfram x- och y-värden som fkn av kurvlängd
S=linspace(0, 1 ,1000);
X=polyval(PX_mean,S);
Y=polyval(PY_mean,S);
%plotta
plot(X-mean(X),Y-mean(Y), 'K', 'linewidth',10)
axis equal

%%

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';


save(filnamn{fil}, 'px', 'py', 'PX_mean', 'PY_mean', '-mat')
disp('save successfull')



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
































