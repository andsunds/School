%% Undersök normalen * avvikelsen, vektoriserat
clc;clf;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom'; 
filnamn{2}='confined_32min_polynom';
filnamn{3}='nonconfined_5min_polynom';
filnamn{4}='nonconfined_167min_polynom';

fil=4;
load(['data/', filnamn{fil}, '.mat'])

N=size(px, 1);

t=(0:(N-1)).'/framerate;
%Bara massa punkter för att få en jämn sträng att hitta TP och mäta längd
S=linspace(0,1,1000);

n=100;%antalet punkter att kolla korr. i

Kn=zeros(N,1);%init

tic
P=linspace(0,1,n);%loopa över hela strängen

%Beräkna normalvektorer i alla punkter och tider. 
[~, N0]=tangent_normal(PX_mean, PY_mean, P);

Q0=[polyval(PX_mean, P); polyval(PY_mean, P)];%Den undersökt punkten på medelkurvan
Q0=bsxfun(@minus, Q0, mean([polyval(PX_mean, S); polyval(PY_mean, S)], 2) );%tyngdpunkt i origo

for i=1:N %loopa över alla bilder (tid)
    Q1=[polyval(px(i,:), P); polyval(py(i,:), P)]; 
    Q1=bsxfun(@minus, Q1, mean([polyval(px(i,:), S); polyval(py(i,:), S)], 2) );%tyngdpunkt i origo
    
    %       v-- specialare som kollar alla normaler
    Kn(i)  = sum(diag(N0.'*(Q1-Q0))); %kolla normalen
end
toc

Kn=Kn/n;%medelvärde av summan
%Strängens totala längd:
L=arclength(PX_mean, PY_mean);

plot(t,Kn)

typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel

xlabel('tid /[s]', 'interpreter', 'LaTeX', 'fontsize',25)
ylabel('$(\mathbf{r}-\mathbf{r}_0)\cdot\mathbf{\hat{n}}_0$ /[px]', 'interpreter', 'LaTeX', 'fontsize',25)
grid on


%F=fft(Kn);
%plot(abs(F(1:end/2)))

%% Jämför tangentvektorerna momentant med medelsträngen
clc;clf;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom'; 
filnamn{2}='confined_32min_polynom';
filnamn{3}='nonconfined_5min_polynom';
filnamn{4}='nonconfined_167min_polynom';

fil=3;
load(['data/', filnamn{fil}, '.mat'])

N=size(px, 1);
t=(0:(N-1)).'/framerate;
S=linspace(0,1,1000);




Kt=zeros(N,1);%init

n=100;%antalet punkter att kolla korr. i
P=linspace(0,1,n);

%Beräknar tangentvektorer i alla punkter på medelsträngen 
T0=tangent_normal(PX_mean, PY_mean, P);

tic
for i=1:N
    T1=tangent_normal(px(i,:), py(i,:), P);
    
    %                 % v-- Skalärprudukt mellan alla möjliga tangentvektorer
    Kt(i)=sum(diag((T1.'*T0), 1))/n;%men det är bara diagonalerna som är intressanta
    
end
toc

plot(t,Kt)
grid on
typ=regexp(filnamn{fil}, '_\d+', 'split');%plockar ut strängtypen
title(sprintf('Fil nr: %d (%s)', fil, typ{1}))%titel
xlabel('tid /[s]', 'interpreter', 'LaTeX', 'fontsize',25)
ylabel('$\mathbf{\hat{t}}_1\cdot\mathbf{\hat{t}}_0$ ', 'interpreter', 'LaTeX', 'fontsize',25)








































