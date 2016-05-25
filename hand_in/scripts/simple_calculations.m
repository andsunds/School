%%
clc;
%clf;
clearvars

load('constants.mat')

t_halv=16.2*3600; %s
lambda=log(2)/t_halv; %s^(-1)

t_0=54*3600; %s
A_dos=2.1e6; %s^(-1)
i=10e-6; %A

sigma=33e-3*1e-24;%M^2
M=76;%g/mol

N_0=5*A_dos/lambda*exp(lambda*t_0) % stycken

n=5*A_dos*e/(i*sigma)*exp(2*lambda*t_0)/(exp(lambda*t_0)-1) %kärnor/cm^2


rho=M*n/N_A 
