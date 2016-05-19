%%
clc;
%clf;
clearvars

t_halv=16.2*3600; %s
lambda=log(2)/t_halv; %s^(-1)

t_0=54*3600; %s
A_dos=2.1e6; %s^(-1)

N_0=5*A_dos/lambda*exp(lambda*t_0) % 
