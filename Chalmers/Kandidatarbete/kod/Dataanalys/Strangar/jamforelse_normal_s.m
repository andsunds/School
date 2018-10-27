%Plotta upp MSD för normalavståndet och avståndet mellan samma s-koordinat
%samt jämför de båda metoderna

%Använder kod från jamvikt.m

clc,clear all
figure(4),clf,figure(3),clf
for fil=1:4
clearvars -except fil

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
XP=zeros(N,res);
YP=zeros(N,res);
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
%k_mt=(Y_m(Q+1)-Y_m(Q-1))/(X_m(Q+1)-X_m(Q-1));
%k_n=-1/k_mt;
b=Y_m(Q)-k_n*X_m(Q);
Nx=k_n.*x+b;

for i=1:N
dXP=polyder(px(i,:));%derivera polynomet i x
dYP=polyder(py(i,:));%derivera polynomet i y

%Nor=[-polyval(dYP, punkt), polyval(dXP, punkt) ]; %normalvektor
Tan=[ polyval(dXP, punkt), polyval(dYP, punkt) ]; %tangentvektor
%Nor=Nor/norm(Nor);  %normerad normalvektor till de olika strängarna
Tan=Tan/norm(Tan);  %normerad tangentvektor

%Tx = tangent till jmvkt-poly som funktion i x-led Tx=k_t*x+m
k_t=Tan(2)/Tan(1);
m=YP(i,Q)-k_t*XP(i,Q);
Tx=k_t.*x+m;

%Approximation till att finna där normalen skär strängen
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

avs(i,1)=XP(i,D(i))-X_m(Q);   %avstånd från jämviktsläge i x-led
avs(i,2)=YP(i,D(i))-Y_m(Q);   %avstånd från jämviktsläge i y-led
end



%MSD, från avs ovan mätt från en punkt=Q samt normalt mot jämvikt
%För att jämföra metoden med att välja samma s eller gå normalt ut från
%jämviktssträngen

MSD_n=zeros(N,2);
MSD_Q=zeros(N,2);
figure(3)
subplot(2,2,fil)
MSD_n=avs(:,1).^2+avs(:,2).^2; %MSD normalt ut från jämviktssträngen
plot(1:N,MSD_n)
figtitle('MSD för en punkt') %Funktion i separat fil
MSD_Q=(XP(:,Q)-X_m(Q)).^2+(YP(:,Q)-Y_m(Q)).^2; %MSD för samma s-koordinat
hold on
plot(1:N,MSD_Q)
hold off
legend('Normalt avstånd', 'För samma s', 'Location','Best')

figure(4)
subplot(2,2,fil);
Diff_t_Q=sqrt(MSD_n)-sqrt(MSD_Q); %Skillnaden mellan dem, ej MSD
plot(1:N,Diff_t_Q)
figtitle('Diff normalt mot samma s')



end