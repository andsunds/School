%% s(dt) = 1/(#particles) * sum( ((x(t)-x(0)).^2 ) over all particles
clc; clearvars
figure(2);clf;pause(.1)

load('../kompleterande_data.mat',...
     'filnamn', 'intensitet', 'std_n', 'std_t', 'sigma_brus', '-mat')

for fil=1:2;
data =load(['../', filnamn{fil}]);
C = separera(data);

%kortare namn
I=intensitet{fil};


%Sammanvägd koef för T och N ger rörligheter som går att jämföra
lambda=sqrt(std_n{fil}.^2+std_t{fil}.^2 - 2*sigma_brus(fil)^2 );
%beräkna koefficienter för lutning i loglog
koef=[ones(size(I)), log(I)]\log(lambda);
koef(1)=exp(koef(1));%konverterar till potenssamband


%initialisering
s=zeros(1000,2);

tic
%loop över alla partiklar med tidsdata upp till 10s
index=find(cellfun('length',C)==1000).'; 
for i=index;
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    tmp=(TN.^2);
    
    %bygger "medelvärde"
    s=s+tmp/(koef(1)*I(i).^koef(2));%normerat medelvärde
    
end
toc

s=sum(s,2)/length(index);


start=10;
stop=101;
t=C{i}(:,1);
c=[ones(stop+1-start,1), log(t(start:stop))]\log(s(start:stop,:));

x=logspace( -2, log10(t(end))).';
y=repmat(exp(c(1,:)), length(x),1).* bsxfun(@power, x, c(2,:));


%plottar data
subplot(1,2,fil)
plot(t, s), hold on
plot(x,y)

str1=sprintf('%.1d t^{%1.2f}', exp(c(1,1)), c(2,1));
%str2=sprintf('%.1d t^{%1.2f}', exp(c(1,2)), c(2,2));


%legend('T', 'N', str1, str2, 'location', 'NorthWest')
legend('Data', str1, 'location', 'NorthWest')
title(filnamn{fil}(1:end-4))
xlabel('Tid/[s]', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15,'XScale','log','YScale','log');
pause(.1)
end

%% S(dt)=(1/T) sum((f(t)-f(t+dt)).^2) over all t
clc; clearvars
figure(3);clf;
%figure(4);clf;
pause(.1)

load('../kompleterande_data.mat',...
     'filnamn', 'intensitet', 'std_n', 'std_t', 'sigma_brus', '-mat')

N=1000;
DT=(1:N).';
addpath('../../');%Lägger till så att create_indecis kan användas
INDECIES=create_indecis(N);
LENGHTS=N+1-DT;


for fil=1:2;
data =load(['../', filnamn{fil}]);
C = separera(data);

%kortare namn
I=intensitet{fil};

%Detta är just nu samma normering som för rörligheten...
lambda=sqrt(std_n{fil}.^2+std_t{fil}.^2 -2*sigma_brus(fil).^2 );
%beräkna koefficienter för lutning i loglog
koef=[ones(size(I)), log(I)]\log(lambda);
koef(1)=exp(koef(1));%konverterar till potenssamband

S=zeros(N,2);
tic
index=find(cellfun('length',C)==N).'; 
for i=index
    TN=koordinatbyte(C{i}(:,2:3));%laddar in data för partikeln
    
    tmp=zeros(length(DT),2);
    
    %Hög minnesåtgång
    A=repmat(TN(:,1), 1,N);
    A=triu(A.'-A);
    tmp(:,1)=sum(A(INDECIES).^2, 2)./LENGHTS;
    
    A=repmat(TN(:,2), 1,N);
    A=triu(A.'-A);
    tmp(:,2)=sum(A(INDECIES).^2, 2)./LENGHTS;
    
    
    % Detta är samma normering som för rörligheten, 
    % kanske skulle man hitta på något annat.
    S=S+(tmp)/(koef(1)*I(i).^koef(2));
end
toc

S=sum(S,2)/length(index);

start=10; show=101; %hur många tidssteg ska undersökas

%anpassar exponentialsamband
Dt=(DT-1)*1e-2;%verklig tid

c=[ones(show+1-start,1), log(Dt(start:show))]\log(S(start:show,:));

x=logspace(-2, log10(Dt(end)) ).';
y=repmat(exp(c(1,:)), length(x),1).* bsxfun(@power, x, c(2,:));



% plottar data
figure(3);
subplot(1,2,fil)
plot(Dt,S), hold on
plot(x,y)


%axis([0, Dt(end), 0, 6e-6]);

str1=sprintf('%.1d dt^{%1.2f}', exp(c(1,1)), c(2,1));
%str2=sprintf('%.1d dt^{%1.2f}', exp(c(1,2)), c(2,2));

%legend('T', 'N', str1,str2, 'location', 'NorthWest')%Om man seprarerar T och N
legend('T', str1, 'location', 'NorthWest')
title(filnamn{fil}(1:end-4))
xlabel('Tidssteg dt/[s]', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('S(dt)', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15,'XScale','log','YScale','log');
pause(.1)

end


%% Båda samtidigt
clc;clearvars
figure(4), pause(.1);clf

load('../kompleterande_data.mat', 'filnamn', '-mat')

N=1000;
Dt=(0:(N-1)).'*1e-2;

data_S=zeros(N,3);
data_c=zeros(2,2);

for fil=1:2;
    tic
    s=MSD_s(fil, N);
    S=MSD_S(fil, N);
    toc
    
    
    %Sparar data för sep. plottning
    data_S(:,fil+1)=S;


    
    subplot(1,2,fil)
    plot(Dt,S), hold on
    plot(Dt,s)

    l=legend('$S$','$s$');
    set(l, 'Interpreter', 'Latex','FontSize',15)

    title(filnamn{fil}(1:end-4))
    xlabel('Tidssteg dt/[s]', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
    ylabel('MSD', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
    set(gca, 'FontSize' ,15,'XScale','log','YScale','log');
    pause(.1)

end

%data_S(:,1)=Dt;

%save('MSD_S.tsv', 'data_S', '-ascii', '-tabs')
%save('MSD_S_linreg.tsv', 'data_S', '-ascii', '-tabs')

