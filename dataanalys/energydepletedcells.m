%% Energydepletedcells

X = csvread('energydepletedcells.csv');

[C index] = separera2(X);




%% 
B= [1;X(:,4)];
P = [X(:,4);1];
E = B-P;
find(E)
%% Testa att C{i} har samma intensitet
 
for i = 1:302
    S = size(C{i});
    x = linspace(0,1,S(1));
    plot(x,C{i}(:,4))
    hold on
end


