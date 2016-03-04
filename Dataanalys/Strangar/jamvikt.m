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
