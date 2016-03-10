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
S=linspace(0, 1 ,res); %Evenly spaced S-coordinates from 0 to 1
X_m=polyval(PX_mean,S); %Evaluate mean x-polynomial at S
Y_m=polyval(PY_mean,S);
X_m=X_m-mean(X_m); %To place center at the origin
Y_m=Y_m-mean(Y_m);
N=size(px,1);
for i=1:N
    XP(i,:)=polyval(px(i,:),S);
    YP(i,:)=polyval(py(i,:),S);
    XP(i,:)=XP(i,:)-mean(XP(i,:));
    YP(i,:)=YP(i,:)-mean(YP(i,:));
end
punkt=.5; %punkt vi kollar på i intervallet [0,1]
Q=length(S)*punkt; %Index of chosen s-cordinate
avs(:,1)=XP(:,Q)-X_m(Q); %x distance to mean string with the same s-coordinate
avs(:,2)=YP(:,Q)-Y_m(Q);
avs=avs';

%Tangenten till jmvkt läge
dX_m=polyder(PX_mean);%derivera jämviktspolynomet i x
dY_m=polyder(PY_mean);%derivera jämviktspolynomet i y

jmN=[-polyval(dY_m, punkt), polyval(dX_m, punkt) ]; %normalvektor till jmvkt
jmT=[ polyval(dX_m, punkt), polyval(dY_m, punkt) ]; %normalvektor till jmvkt
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
X_m=polyval(PX_mean,S);
Y_m=polyval(PY_mean,S);
X_m=X_m-mean(X_m);
Y_m=Y_m-mean(Y_m);
N=size(px,1);
x=linspace(min(X_m),max(X_m),res);
for i=1:N
    XP(i,:)=polyval(px(i,:),S);
    YP(i,:)=polyval(py(i,:),S);
    XP(i,:)=XP(i,:)-mean(XP(i,:));
    YP(i,:)=YP(i,:)-mean(YP(i,:));
end
punkt=.5; %punkt vi kollar på i intervallet [0,1]
Q=length(S)*punkt;

%Tangenten till jmvkt läge
dX_m=polyder(PX_mean);%derivera jämviktspolynomet i x
dY_m=polyder(PY_mean);%derivera jämviktspolynomet i y

jmN=[-polyval(dY_m, punkt), polyval(dX_m, punkt) ]; %normalvektor till jmvkt
jmT=[ polyval(dX_m, punkt), polyval(dY_m, punkt) ]; %tangentvektor till jmvkt
jmN=jmN/norm(jmN);  %normerad normalvektor
jmT=jmT/norm(jmT);  %normerad tangentvektor

%Nx = normal till jmvkt-poly som funktion i x-led Nx=k_n*x+b
k_n=jmN(2)/jmN(1);
b=Y_m(Q)-k_n*X_m(Q);
Nx=k_n.*x+b;

for i=1:N
dXP=polyder(px(i,:));%derivera polynomet i x
dYP=polyder(py(i,:));%derivera polynomet i y

Nor=[-polyval(dYP, punkt), polyval(dXP, punkt) ]; %normalvektor
Tan=[ polyval(dXP, punkt), polyval(dYP, punkt) ]; %tangentvektor
Nor=Nor/norm(Nor);  %normerad normalvektor
Tan=Tan/norm(Tan);  %normerad tangentvektor

%Tx = tangent till jmvkt-poly som funktion i x-led Tx=k_t*x+m
k_t=Tan(2)/Tan(1);
m=YP(i,Q)-k_t*XP(i,Q);
Tx=k_t.*x+m;

C=find(min(abs(Nx-Tx))==abs(Nx-Tx));  %skärning mellan tangent och normal
val=C/res*(max(X_m)-min(X_m))+min(X_m);  %hittar x-värdet på polynomet i skärningspunkten
D(i)=find(min(abs(XP(i,:)-val))==abs(XP(i,:)-val)); %skärningspunkten på polynomet


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

%Plotta ut beräknade punkter och normal/tangent
% plot(x,Nx,'r',x,Tx,'r') %Normal and tangent
% hold on
% plot(X_m,Y_m, XP(i,:),YP(i,:)) %Mean string and current string
% plot(X_m(Q),Y_m(Q),'o',XP(i,Q),YP(i,Q),'*',XP(i,D),YP(i,D),'o') %Chosen s-coordinate on mean and current string and intersection point
% axis([min(x)-10 max(x)+10 min(Y_m)-20 max(Y_m)+20])
% plot(x,abs(Tx-Nx))
% hold off
% pause(0.3)

avs(i,1)=XP(i,D(i))-X_m(Q);   %avstånd från jämviktsläge i x-led
avs(i,2)=YP(i,D(i))-Y_m(Q);   %avstånd från jämviktsläge i y-led
end

L=(sum(sqrt(sum(diff([polyval(PX_mean, S); polyval(PY_mean, S)] ,2).^2 ,1)))); %längd på sträng
Nprod=jmN*avs'/L;             %skalärprod mellan avstånd och normalvektor, normerad med stärngens längd.
int=zeros(1,N);
 
 
for k=1:N
    int(k)=sum(Nprod(1:k))+Nprod(k);   %summerar avvikelser fram till punkten k, motsvarar integral.
end
 
 
int=int/N;
figure(2)
subplot(1,2,1),plot(1:N,Nprod,1:N,0*[1:N],'-'),title('Avvikelse från jämviktsläge')
subplot(1,2,2),plot(int),title('summa av avvikelser (normerad med N)')


end
