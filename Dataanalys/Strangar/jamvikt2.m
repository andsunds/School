%% Undersök normalen * avvikelsen
clc;clf;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=1;
load(['data/', filnamn{fil}])

N=size(px, 1);
S=linspace(0,1,1000);

%Strängens totala längd
L=sum(sqrt(sum(diff([polyval(PX_mean, S); polyval(PY_mean, S)] ,2).^2 ,1)))

n=100;%antalet punkter att kolla korr. i

Kn=zeros(N,1);%init
%Kt=zeros(N,1);%init

for P=linspace(0,1,50)%loopa över hela strängen

%P=0.5;
eps=0;%Kolla medelvärde av kurvsegment
p=P+linspace(-eps,eps,n);%Vilka punkter vi ska kolla efter tangentvektor

%Beräkna tangentvektorer i alla punkter och tider. 
[T0, N0]=tangent_normal(PX_mean, PY_mean, P);

Q0=[polyval(PX_mean, P); polyval(PY_mean, P)];%Den undersökt punkten på medelkurvan
Q0=Q0-mean([polyval(PX_mean, S); polyval(PY_mean, S)], 2);%tyngdpunkt i origo


for i=1:N
    Q1=mean([polyval(px(i,:), p); polyval(py(i,:), p)], 2); %medelvärdet över kurvsegmentet
    Q1=Q1-mean([polyval(px(i,:), S); polyval(py(i,:), S)], 2);%tyngdpunkt i origo
    
    %Kn(i)=Kn(i)+N0.'*(Q1-Q0);
    Kn(i)=Kn(i)+sqrt((Q1-Q0).'*(Q1-Q0));
end
end

plot(Kn/50)

%F=fft(Kn);
%plot(abs(F(1:end/2)))

%% Jämför tangentvektorerna momentant med medelsträngen
clc;clf;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=3;
load(['data/', filnamn{fil}])

N=size(px, 1);
S=linspace(0,1,1000);




Kt=zeros(N,1);%init

n=100;%antalet punkter att kolla korr. i
P=linspace(0,1,n);

%Beräknar tangentvektorer i alla punkter på medelsträngen 
T0=tangent_normal(PX_mean, PY_mean, P);

tic
for i=1:N
    T1=tangent_normal(px(i,:), py(i,:), P);
    
    %             % v-- Skalärprudukt mellan alla möjliga tangentvektorer
    Kt(i)=mean(diag((T1.'*T0), 1));%men det är bara diagonalerna som är intressanta
    
end
toc

plot(Kt)









































