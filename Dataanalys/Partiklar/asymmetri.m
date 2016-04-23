%Asymmetri

%Simulering vanlig Brownsk, sannolikhet för olika asymmetrier då man mäter på 193
%partiklar
clc;clf; clearvars
hold on
N_steps=1000; %Antal steg för varje partikel
N_trials=193; %antal partiklar
M_trials=300; %antal mätserier med N_trials antal partiklar

Asym_mean_sim=zeros(M_trials,1);

for j=1:M_trials
    
Asym_sim_upper=zeros(N_trials,1); %För att ta medelvärdet av senare i täljaren
Asym_sim_lower=zeros(N_trials,1); %För att ta medelvärdet av senare i nämnaren
    
for i=1:N_trials
    X=[0; cumsum(randn(N_steps-1,1))]; %Normalfördelade steg, start i 0
    Y=[0; cumsum(randn(N_steps-1,1))];
    gyr = [sum(X.^2)/N_steps-(sum(X)/N_steps)^2, sum(X.*Y)/N_steps-(sum(Y)/N_steps)*(sum(X)/N_steps); 
           sum(Y.*X)/N_steps-(sum(Y)/N_steps)*(sum(X)/N_steps), sum(Y.^2)/N_steps-(sum(Y)/N_steps).^2];
    [~,D] = eig(gyr); %Egenvrden läng diagonalen: D(1) och D(4)

    Asym_sim_upper(i)=(D(1)-D(4))^2;
    Asym_sim_lower(i)=(D(1)+D(4))^2;
    %Asym_sim_upper(i)=(D(1)^2+D(4)^2-D(1)*D(4)); %Alternativ definition
    %Asym_sim_lower(i)=(D(1)+D(4))^2;
end

Asym_mean_sim(j)=mean(Asym_sim_upper)/mean(Asym_sim_lower); %A för aktuell simulering

end

hist(Asym_mean_sim) %Fördelingen för A
%median(Asym_mean_sim)
mean(Asym_mean_sim) %Att jämföra med tabellvärde nedan

d=2; %dimension
Asym_litterature=2*(d+2)/(5*d+4) %Då N_trials går mot oändligheten

%%
%Asymmetrin för datan

load('filnamn.mat')


Asym_data_mean=[0,0]; %Medelvärde för energydepleted och logphase

for fil=1:2 %För att undersöka datan för båda celltyperna
    
%laddar in data
data =load(filnamn{fil});
C = separera(data);
n=length(C);%antal partiklar

Asym_data_upper=zeros(n,1); %För att ta medelvärdet av senare i täljaren
Asym_data_lower=zeros(n,1); %För att ta medelvärdet av senare i nämnaren

for i=1:n
    X=C{i}(:,2);
    Y=C{i}(:,3);
    N_steps=length(X);
    gyr = [sum(X.^2)/N_steps-(sum(X)/N_steps)^2, sum(X.*Y)/N_steps-(sum(Y)/N_steps)*(sum(X)/N_steps); 
           sum(Y.*X)/N_steps-(sum(Y)/N_steps)*(sum(X)/N_steps), sum(Y.^2)/N_steps-(sum(Y)/N_steps).^2];
    [~,D]=eig(gyr);
    
    Asym_data_upper(i)=(D(1)-D(4))^2;
    Asym_data_lower(i)=(D(1)+D(4))^2;
    %Asym_data_upper(i)=(D(1)^2+D(4)^2-D(1)*D(4)); %Alternativ definition
    %Asym_data_lower(i)=(D(1)+D(4))^2;
end

Asym_data_mean(fil)=mean(Asym_data_upper)/mean(Asym_data_lower) %Medelvärde för energydepleted och logphase

end