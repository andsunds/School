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

%% 2-korr (kör övre koden först)
clearvars -except Pos AccPos prod fil
fil=1;
avs=prod{fil,1}+prod{fil,2};
N=length(avs);
M=floor(N/2);

s = 80;
t = 20;
normering=(s+1)*(t+1);
G = zeros(2*s+1,t+1);

for z = M-floor(s/2):M+ceil(s/2)
    for tau = 0:t
        for k = -s:s
            for j = 0:t
            G(z+floor(s/2)+1,tau+1) = G(z+floor(s/2)+1,tau+1)+...
                +XY{j+1}(z,2)*XY{tau+j+1}(z,2)/normering;
            end
        end
    end
end
       


%% plottar

clf
i=1;
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
%% avvikelse definerad som skärning mellan jämviktssträngens tangentlinje och sträng
clc,clear all
for fil=4:4

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
x=linspace(min(X),max(X),res);
for i=1:N
    XP(i,:)=polyval(px(i,:),S);
    YP(i,:)=polyval(py(i,:),S);
    XP(i,:)=XP(i,:)-mean(XP(i,:));
    YP(i,:)=YP(i,:)-mean(YP(i,:));
end
punkt=.5; %punkt vi kollar på i intervallet [0,1]
Q=length(S)*punkt;

%Tangenten till jmvkt läge
dX=polyder(PX_mean);%derivera jämviktspolynomet i x
dY=polyder(PY_mean);%derivera jämviktspolynomet i y

jmN=[-polyval(dY, punkt), polyval(dX, punkt) ]; %normalvektor till jmvkt
jmT=[ polyval(dX, punkt), polyval(dY, punkt) ]; %tangentvektor till jmvkt
jmN=jmN/norm(jmN);  %normerad normalvektor
jmT=jmT/norm(jmT);  %normerad tangentvektor

b=Y(Q)-jmN(2)/jmN(1)*X(Q);
Nx=jmN(2)/jmN(1).*x+b; %normal till jmvkt-poly i som funktion i x-led

for i=1:N
dXP=polyder(px(i,:));%derivera polynomet i x
dYP=polyder(py(i,:));%derivera polynomet i y

Nor=[-polyval(dYP, punkt), polyval(dXP, punkt) ]; %normalvektor
Tan=[ polyval(dXP, punkt), polyval(dYP, punkt) ]; %tangentvektor
Nor=Nor/norm(Nor);  %normerad normalvektor
Tan=Tan/norm(Tan);  %normerad tangentvektor


m=YP(i,Q)-Tan(2)/Tan(1)*XP(i,Q);
Tx=Tan(2)/Tan(1).*x+m;

C=find(min(abs(Nx-Tx))==abs(Nx-Tx));                    %skärning mellan tangent och normal
val=C/res*(max(X)-min(X))+min(X);                       %hittar x-värdet på polynomet i skärningspunkten
D(i)=find(min(abs(XP(i,:)-val))==abs(XP(i,:)-val));     %skärningspunkten på polynomet


%liten loop som kastar orimliga punkten på polynomet. Gäller för strängar
%som har 2 värden för ett x-värde. Dvs kröker runt sig själv i x-led.
if i>1 && D(i)-D(i-1)>res/10
    for e=1:1000
    E=find(abs(XP(i,:)-val)<1/(e^2));
    if length(E)<10
        F=find(min(abs(E-Q))==abs(E-Q));
        D(i)=E(F);
        break
    end
    end
end

%plot(x,Nx,X,Y,XP(i,:),YP(i,:),X(Q),Y(Q),'o',XP(i,Q),YP(i,Q),'*',x,Tx,XP(i,D(i)),YP(i,D(i)),'o')
%axis equal
%plot(x,abs(Tx-Nx))
avs(i,1)=XP(i,D(i))-X(Q);       %avstånd från jämviktsläge i x-led
avs(i,2)=YP(i,D(i))-Y(Q);       %avstånd från jämviktsläge i y-led
end

L=(sum(sqrt(sum(diff([polyval(PX_mean, S); polyval(PY_mean, S)] ,2).^2 ,1)))); %längd på sträng
Nprod=jmN*avs'/L;                %skalärprod mellan avstånd och normalvektor, normerad med stärngens längd.
int=zeros(1,N);
 
 
for k=1:N
    int(k)=sum(Nprod(1:k))+Nprod(k);   %summerar avvikelser fram till punkten k, motsvarar integral.
end
 
 
int=int/N;
figure(2)
subplot(1,2,1),plot(1:N,Nprod,1:N,0*[1:N],'-'),title('Avvikelse från jämviktsläge')
subplot(1,2,2),plot(int),title('summa av avvikelser (normerad med N)')


end












