%% Korrelation i rörelsen
clc;clf;clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

offs=[0.007, 0.009];

for fil=1:2 %för att plotta datan för båda celltyperna

data =load(filnamn{fil});


C = separera(data);
n=length(C);%antal partiklar

storl=zeros(n,1);%stoleken på partiklarna
medelsteg=zeros(n,1);%medelstegen på en prtikel
steg_osakerhet=zeros(n,1);%avvikelen från medelsteglängden
rorlighet=zeros(n,1);%hur stort område vandrar partikeln runt i

hold on

for i=1:n
    %if size(C{i},1)~=1000
    %    continue
    %end
    XY=C{i}(2:end-1,2:3);
    D1=diff(XY(1:end-1));
    D2=diff(XY(2:end));
    D=diff(XY);
    
    pause(0.01)
    subplot(2,2,2*fil-1)
    plot(XY(1:end-1,:),D, '.r');hold on;axis equal
    
    subplot(2,2,2*fil)
    plot(D1,D2, '.b');hold on;axis equal
end




end