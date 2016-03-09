%% Vektoriserad tangentvektors korrelation
%<t(s) * t(s+l)> ~ exp(-l/L_P)?
% s är en parameter som går från 0 till 1 oberoende av strängens längd.
clc;clf;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom'; 
filnamn{2}='confined_32min_polynom';
filnamn{3}='nonconfined_5min_polynom';
filnamn{4}='nonconfined_167min_polynom';

fil=3;
load(['data/', filnamn{fil}, '.mat'])

N=size(px, 1);

n=100;%antalet punkter att kolla korr. i
K=zeros(n,1);%init.
s=linspace(0,1,n);%Vilka punkter vi ska kolla efter tangentvektor

%Beräkna tangentvektorer i alla punkter och tider. 
tangent=tangent_normal(px, py, s);

addpath('../');%Lägger till så att create_indecis kan användas
INDEX=create_indecis(n);%Tar fram index som sorterar längs med diagonalerna

tic
for i=1:N; %loopa över alla bilder
%                        /-- tar fram tangentvektorerna vid just tid i och
%   /-- bara övre        |    matrismultiplicerar dem för att få en matris
%   V    diagonalerna    V      med alla produkter
tmp=triu(tangent(:,:,i).'*tangent(:,:,i));

%  v-- lägger bidrag från varje bild  %v-- medelvärde över alla bilder
K=K+sum( tmp(INDEX), 2);
%                       ^normerar m.a.p. antalet datapunkter
end
K=K./fliplr(1:n).'/N;%medelvärde i tid
toc

% Plottning 

%Strängens totala längd:
S=linspace(0,1,1000);
L=sum(sqrt(sum(diff([polyval(PX_mean, S); polyval(PY_mean, S)] ,2).^2 ,1)));
l=L*s;%omvandla från paramaterna s, till längden på strängen
plot(l, K), hold on


% Anpassa exponentialkurva
i_bra=(floor(n*0.05 +1):floor(n*.5))';

C=[l(i_bra)' ones(size(i_bra))]\log(K(i_bra));

L_P=-1/C(1);
K_0=exp(C(2));

plot(l, K_0*exp(-l/L_P))

str=sprintf('$%.3f \\exp(-\\Delta l/%3.f)$', K_0, L_P);
leg=legend('Ber\"a{}knad korrelation', str);
set(leg, 'interpreter', 'Latex')

xlabel('$\Delta l$, l\"a{}ngs str\"a{}ngen','Interpreter','Latex');
ylabel('$<\mathbf{t}(l)\cdot\mathbf{t}(l+\Delta l)>_{l, t}$','Interpreter','Latex')
set(gca,'Fontsize',24)%, 'yscale','log');




%% Tangentvektorskorr i tid 
clf;clc;clear all

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom'; 
filnamn{2}='confined_32min_polynom';
filnamn{3}='nonconfined_5min_polynom';
filnamn{4}='nonconfined_167min_polynom';

fil=2;
load(['data/', filnamn{fil}, '.mat'])

N=size(px, 1); % Total tid

n=300;%antalet punkter att kolla korr. i
K=zeros(N,n);%init.
eps=.1;
s=linspace(eps,1-eps,n);%Vilka punkter vi ska kolla efter tangentvektor

%Beräkna tangentvektorer i alla punkter och tider. 
tangent=tangent_normal(px, py, s);

for i=1:N
    T1 = tangent(:,:,i); % Tangentvektorer tid: i
    for dt=0:(N-i)
        T2 = tangent(:,:,i+dt); % Tangentvektorer tid: i+dt
        K(dt+1,:) = K(dt+1)+sum(T1.*T2, 1)./(N-dt);
        %                   ^-- summan i skalärprodukten
    end
end

Ks = sum(K,2)'/n; % Summera �ver alla punkter p� str�ngen
dt = 0:(N-1);

%figure(1)
plot(dt,(Ks))
xlabel('$d\tau$ [tid]','Interpreter','Latex');
ylabel('$<\mathbf{t}(\tau)\cdot\mathbf{t}(\tau+d\tau)>$','Interpreter','Latex')
set(gca,'Fontsize',24);%'xscale','log','yscale','log');

%%
figure(2)
surf(K(:,:)) % Summerar ej �ver punkterna p� str�ngen
xlabel('B�gl�ngd s [l�ngd]');
ylabel('dt [tid]');
zlabel('<t(s)*t(s+dt)>')
% figure(1)






