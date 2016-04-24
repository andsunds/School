%Asymmetri

%Simulering vanlig Brownsk, sannolikhet f??r olika asymmetrier d?? man m??ter p?? 193
%partiklar
clc;clf; clearvars
hold on
N_steps=1000; %Antal steg f??r varje partikel
N_trials=193; %antal partiklar
M_trials=300; %antal m??tserier med N_trials antal partiklar

Asym_mean_sim=zeros(M_trials,1);

for j=1:M_trials
    
Asym_sim_upper=zeros(N_trials,1); %F??r att ta medelv??rdet av senare i t??ljaren
Asym_sim_lower=zeros(N_trials,1); %F??r att ta medelv??rdet av senare i n??mnaren
    
for i=1:N_trials
    X=[0; cumsum(randn(N_steps-1,1))]; %Normalf??rdelade steg, start i 0
    Y=[0; cumsum(randn(N_steps-1,1))];
    gyr = [sum(X.^2)/N_steps-(sum(X)/N_steps)^2, sum(X.*Y)/N_steps-(sum(Y)/N_steps)*(sum(X)/N_steps); 
           sum(Y.*X)/N_steps-(sum(Y)/N_steps)*(sum(X)/N_steps), sum(Y.^2)/N_steps-(sum(Y)/N_steps).^2];
    [~,D] = eig(gyr); %Egenvrden l??ng diagonalen: D(1) och D(4)

    Asym_sim_upper(i)=(D(1)-D(4))^2;
    Asym_sim_lower(i)=(D(1)+D(4))^2;
    %Asym_sim_upper(i)=(D(1)^2+D(4)^2-D(1)*D(4)); %Alternativ definition
    %Asym_sim_lower(i)=(D(1)+D(4))^2;
end

Asym_mean_sim(j)=mean(Asym_sim_upper)/mean(Asym_sim_lower); %A f??r aktuell simulering

end

hist(Asym_mean_sim) %F??rdelingen f??r A
%median(Asym_mean_sim)
mean(Asym_mean_sim) %Att j??mf??ra med tabellv??rde nedan

d=2; %dimension
Asym_litterature=2*(d+2)/(5*d+4) %D?? N_trials g??r mot o??ndligheten

%%
%Asymmetrin f??r datan

load('filnamn.mat')


Asym_data_mean=[0,0]; %Medelv??rde f??r energydepleted och logphase

for fil=1:2 %F??r att unders??ka datan f??r b??da celltyperna
    
%laddar in data
data =load(filnamn{fil});
C = separera(data);
n=length(C);%antal partiklar

Asym_data_upper=zeros(n,1); %F??r att ta medelv??rdet av senare i t??ljaren
Asym_data_lower=zeros(n,1); %F??r att ta medelv??rdet av senare i n??mnaren

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

Asym_data_mean(fil)=mean(Asym_data_upper)/mean(Asym_data_lower) %Medelv??rde f??r energydepleted och logphase

end

%%
%fBm

Asym_fBm=@(H)2-(1/(2*(H+1)^2))/(1/(2*(H+1)^2)+(2*H+1)/(4*(4*H+1))-1/(4*H+3)-(gamma(2*H+2))^2/(gamma(4*H+4)));
Asym_fBm(0.5) %f?r Wienerproces
A_find=@(H,A)abs(Asym_fBm(H)-A); %Skillnad mellan ber?knad och teoretisk

H_data=fminbnd(@(H)A_find(H,Asym_data_mean(1)),0,1)