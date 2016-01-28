%% Energydepletedcells

Y = csvread('logphasecells.csv');

[D index2] = separera(Y);



%% Testa att D{i} har samma intensitet
for i = 1:193
    S = size(D{i});
    x = linspace(0,1,S(1));
    plot(x,D{i}(:,4))
    hold on
end


