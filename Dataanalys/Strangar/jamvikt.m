%% avvikelse från jämvikt m.a.p tangent/normal - vektorer
clc,clear all
AccPos=cell(4,2);
prod=cell(4,2);
for fil=1:4
    clearvars -except Pos AccPos prod fil 
%fil=1; % fil att kolla på

if fil==1
load('data/confined_28min_polynom.mat', '-mat');
elseif fil==2
load('data/confined_32min_polynom.mat', '-mat');
elseif fil==3
load('data/nonconfined_5min_polynom.mat', '-mat');
elseif fil==4
load('data/nonconfined_167min_polynom.mat', '-mat');
end

res=1000; %bestämmer antal punkter på "bågen"
S=linspace(0, 1 ,res);
X=polyval(PX_mean,S);
Y=polyval(PY_mean,S);
X=X-mean(X);
Y=Y-mean(Y);
N=size(px,1);
for i=1:N
    XP(i,:)=polyval(px(i,:),S);
    YP(i,:)=polyval(py(i,:),S);
    XP(i,:)=XP(i,:)-mean(XP(i,:));
    YP(i,:)=YP(i,:)-mean(YP(i,:));
end
punkt=.5; %punkt vi kollar på i intervallet [0,1]
Q=length(S)*punkt;
avs(:,1)=XP(:,Q)-X(Q);
avs(:,2)=YP(:,Q)-Y(Q);
avs=avs';

%Tangenten till jmvkt läge
dX=polyder(PX_mean);%derivera jämviktspolynomet i x
dY=polyder(PY_mean);%derivera jämviktspolynomet i y

jmN=[-polyval(dY, punkt), polyval(dX, punkt) ]; %normalvektor till jmvkt
jmT=[ polyval(dX, punkt), polyval(dY, punkt) ]; %normalvektor till jmvkt
jmN=jmN/norm(jmN);  %normerad normalvektor
jmT=jmT/norm(jmT);  %normerad tangentvektor


for i=1:N;
%avs=avs./repmat(sqrt(sum(avs.^2,1)),2,1);  %normering av avståndsvektorn


%Sprod =skalärprod, N=normalvektor, T=tangentvektor, Acc="integral"
NSprod(i)=jmN*avs(:,i);
TSprod(i)=jmT*avs(:,i);
AccNSprod(i)=sum(NSprod)+jmN*avs(:,i);
AccTSprod(i)=sum(TSprod)+jmT*avs(:,i);

end
prod{fil,1}=TSprod; %skalarprod med tangenten
prod{fil,2}=NSprod; %skalarprod med normalen
AccPos{fil,1}=AccNSprod;
AccPos{fil,2}=AccTSprod;
Pos(fil,1)=sum(NSprod)/sum(abs(NSprod)); %>0 betyder att strängen är mer på positiv sida av jmkvt
Pos(fil,2)=sum(TSprod)/sum(abs(TSprod));
end
%% plottar
clf
i=4;
figure(1),  hold on
subplot(2,3,1),plot(prod{i,1}), title('skalarprod m.a.p normalvektor')
subplot(2,3,2),plot(prod{i,2}), title('skalarprod m.a.p Tangentvektor')
subplot(2,3,3),plot(prod{i,1}+prod{i,2}), title('skalarprod normalvektor+ skalarprod Tangentvektor')
subplot(2,3,4),plot(AccPos{fil,1}), title('summa av skalarprod m.a.p normalvektor')
subplot(2,3,5),plot(AccPos{fil,2}), title('summa av skalarprod m.a.p tangentvektor')
subplot(2,3,6),plot(AccPos{fil,1}+AccPos{fil,2}), title('summa av skalarprod normalvektor+Tangentvektor')

% figure(2), clf
% S=fft(prodT{i});
% plot(abs(S(1:end/2)))
%% avvikelse från jämvikt definerad med skärning med tangentlinje
clc,clear all
for fil=1:1

if fil==1
load('data/confined_28min_polynom.mat', '-mat');
elseif fil==2
load('data/confined_32min_polynom.mat', '-mat');
elseif fil==3
load('data/nonconfined_5min_polynom.mat', '-mat');
elseif fil==4
load('data/nonconfined_167min_polynom.mat', '-mat');
end

res=10000; %bestämmer antal punkter på "bågen"
S=linspace(0, 1 ,res);
X=polyval(PX_mean,S);
Y=polyval(PY_mean,S);
X=X-mean(X);
Y=Y-mean(Y);
N=size(px,1);
for i=1:N
    XP(i,:)=polyval(px(i,:),S);
    YP(i,:)=polyval(py(i,:),S);
    XP(i,:)=XP(i,:)-mean(XP(i,:));
    YP(i,:)=YP(i,:)-mean(YP(i,:));
end
punkt=.99; %punkt vi kollar på i intervallet [0,1]
Q=length(S)*punkt;

%Tangenten till jmvkt läge
dX=polyder(PX_mean);%derivera jämviktspolynomet i x
dY=polyder(PY_mean);%derivera jämviktspolynomet i y

dXP=polyder(XP(1,:));%derivera polynomet i x
dYP=polyder(YP(1,:));%derivera polynomet i y

jmN=[-polyval(dY, punkt), polyval(dX, punkt) ]; %normalvektor till jmvkt
jmT=[ polyval(dX, punkt), polyval(dY, punkt) ]; %tangentvektor till jmvkt
jmN=jmN/norm(jmN);  %normerad normalvektor
jmT=jmT/norm(jmT);  %normerad tangentvektor

Nor=[-polyval(dYP, punkt), polyval(dXP, punkt) ]; %normalvektor
Tan=[ polyval(dXP, punkt), polyval(dYP, punkt) ]; %tangentvektor
Nor=Nor/norm(Nor);  %normerad normalvektor
Tan=Tan/norm(Tan);  %normerad tangentvektor

clf

x=linspace(min(X),max(X),res);
b=Y(Q)-jmN(2)/jmN(1)*X(Q);
Nx=jmN(2)/jmN(1).*x+b;

m=YP(1,Q)-Tan(2)/Tan(1)*XP(1,Q);
Tx=Tan(2)/Tan(1).*x+m;
%C=find(abs(Nx-YP(1,:)-XP(1,:))<=0.2)
plot(x,Nx,x,Y,XP(1,:),YP(1,:),X(Q),Y(Q),'o',XP(1,Q),YP(1,Q),'*',x,Tx)
axis equal

end

















