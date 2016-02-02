%% Steglängd och rörlighet(std) som fkn av storlek
clc;clf;clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

offs=[0.007, 0.009]*0;%det finns lite offest i steglängden

for fil=1:2 %för att plotta datan för båda celltyperna

data =load(filnamn{fil});


C = separera(data);
n=length(C);%antal partiklar

storl=zeros(n,1);%stoleken på partiklarna
medelsteg=zeros(n,1);%medelstegen på en prtikel
steg_osakerhet=zeros(n,1);%avvikelen från medelsteglängden
rorlighet=zeros(n,1);%hur stort område vandrar partikeln runt i

for i=1:n
    storl(i)=mean(C{i}(:,4));
    
    XY=C{i}(:,2:3);
    
    steglangd=sqrt(sum(diff(XY).^2,2));%alla steglängder
    medelsteg(i)=mean(steglangd);
    steg_osakerhet(i)=std(steglangd)/sqrt(length(steglangd));
    
    rorlighet(i)=sqrt(sum(var(XY)));
end

subplot(2,2,2*fil-1)
plot(storl, rorlighet, '.r')
title(filnamn{fil})
xlabel('Storlek', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('R\"o{}rlighet', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15)%,'XScale','log','YScale','log');


subplot(2,2,2*fil)
errorbar(storl, medelsteg-offs(fil), steg_osakerhet, '.')
title(filnamn{fil})
xlabel('Storlek', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
ylabel('Stegl\"a{}ngd', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
set(gca,'FontSize',15)%,'XScale','log','YScale','log');
axis([0,20,0,8e-9])

fit(storl,rorlighet,'power1') %Powerfit r�rlighet
end

%% varians som fkn av tid
clc;clf;clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

offs=[0.007, 0.009];


fil=1;
data =load(filnamn{fil});


C = separera(data);
n=length(C);%antal partiklar

storl=zeros(n,1);%stoleken på partiklarna

VAR_tid=cell(n,1);%hur stort område vandrar partikeln runt i
TID=cell(n,1);

for i=1:n
    storl(i)=mean(C{i}(:,4));
    
    l=length(C{i}(:,1));
    var_tid=zeros(l,1);
    XY=C{i}(:,2:3);
    for j=1:l
        var_tid(j)=(sum(var(XY(1:j))));
    end
    VAR_tid{i}=var_tid;
    TID{i}=C{i}(:,1);
end
fprintf('Långsam for-loop klar\n')
%%
clf
hold on
l=length(TID);
M=zeros(1000,1);
for i=1:l
    pause(0.01);plot(TID{i}, VAR_tid{i}./mean(VAR_tid{i}))
    if length(TID{i})==1000
        M=M+VAR_tid{i}./mean(VAR_tid{i}(1:end));
    end
end

%plot(linspace(0.01,10,1000),M)
%set(gca,'XScale','log','YScale','log');





