%Test för i=59, en av de svagaste intensiteterna, ger ett glapp i y-led på
%ca 2 nm. Mätfel?
clf
clear

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=1;
data =load(filnamn{fil});
C = separera(data);


i=59;
D=C{i};

n=length(D)

%Spela upp partikelns förflyttning
hold on
for i=1:n
    plot(D(i,2),D(i,3),'.g') %Grön prick=aktuell position
    pause(0.2)
    plot(D(i,2),D(i,3),'.r')
end
