%% riktningsberoende
clc;clf;clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';



fil=1;
data =load(filnamn{fil});

C = separera(data);
n=length(C);%antal partiklar

for i=1:n
    v=atan(diff(C{i}(:,3))./(diff(C{i}(:,2))));
    hist(v)
    pause(2)
end




