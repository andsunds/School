%%
clc;clear all;%clf

load Data_GlobalIdx2
prices  = [Dataset.TSX Dataset.CAC Dataset.DAX ...
Dataset.NIK Dataset.FTSE Dataset.SP];

returns =  tick2ret(prices);

nVariables  = size(returns, 2)
expReturn   = mean(returns);
sigma       = std(returns);
correlation = corrcoef(returns);
t           = 0;
X           = 100;
X           = X(ones(nVariables,1));

F = @(t,X) diag(expReturn) * X;
G = @(t,X) diag(X) * diag(sigma);

SDE = sde(F, G, 'Correlation', correlation, 'StartState', X)


[Paths, Times, Z] = simByEuler(SDE, 2)
