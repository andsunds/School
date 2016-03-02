
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
