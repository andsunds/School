%% Some plots 
% Plots E and M as a function of time
clc;clear;clf
beta = 0.1 * (1.096^30);
FILENAME = sprintf( 'bin/beta_%0.5f.bin', beta);
fID=fopen(FILENAME,'rb');
data=fread(fID,[2,inf],'real*8').';
fclose(fID);

figure(1)

subplot(1,2,1)
plot(data(:,1));
title(sprintf('energy E, (beta = %0.5f)',beta));
subplot(1,2,2)
plot(data(:,2));
set(gca, 'yLim',[-1,1])
title(sprintf('order parameter M, (beta = %0.5f)',beta));



%% Extracting data
% Calculates <E> and <M> as a function of temperature, 
% then save that information to  file.
clc;clear;clf

beta0   = 0.25;       % The base value of beta
d_beta  = 16e-3;%1.096;     % The factor the <beta0> is shifted
range   = 0:127;    % The powers of <factor> <beta0> is shifted
t_start = 5e4;       % The index after which we start averaging

NT  = length(range); % # temperature steps
T   = zeros(NT,1);   % init, vector with temperatures
EM  = zeros(NT,2);   % init, matrix with averages of E and M
dEM = zeros(NT,2);   % init, matrix with std in mean of E and M

tic
for i = 1:NT
    % Reading from file:
    beta     = beta0 + (d_beta*range(i));
    FILENAME = sprintf( 'DEBUG/EM_beta_%0.5f.bin', beta);
    fID      = fopen(FILENAME,'rb');
    data     = fread(fID,[2,inf],'real*8').';
    fclose(fID);
    
    % Analyzing the data
    
    T(i)       = 1/beta;
    EM(i,:)    = mean(data(t_start:end, :),    1); 
    stdEM(i,:) =  std(data(t_start:end, :), 0, 1);
    
end
toc 
% This (reading the data and mean + std) took 8 min.
% Generating the data through simulations in C took 30 min.

save_data = [T, EM, stdEM]; %This isn't the same order as from "avg"

%save('T_E_M_dE_dM.tsv', 'save_data', '-ascii', '-double', '-tabs');




%%
clc;clf;clear

data = load('EstdEMstdM_beta_0.400-0.451_2048_PERIODIC.tsv');
%data = load('EstdEMstdM_beta_0.100-2.147_2048_PERIODIC.tsv', '-ascii');
%data = load('EstdEMstdM_beta_0.100-2.147_2048.tsv', '-ascii');
%data = load('TEMstdEstdM_logscale.tsv', '-ascii');




T = data(:,1);
E = data(:,2); stdE = data(:,3)/sqrt(2e7);
M = data(:,4); stdM = data(:,5)/sqrt(2e7);




figure(1)
subplot(1,2,1)
%errorbar(T,E,stdE, 'b.')
plot(T,stdE.^2./T.^2,'.')
%plot(T,E, 'b.')

set(gca, 'xScale', 'lin');
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$\langle E\rangle/J$', 'interpreter', 'LaTeX');



subplot(1,2,2)
errorbar(T,abs(M),stdM, 'b.')
%plot(T,abs(M), 'b.')
plot(T,stdM, '.')

set(gca, 'yScale', 'lin')
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$\langle M\rangle$', 'interpreter', 'LaTeX');







%%
clc;clf;clear

data = load('TEstdErhoXY_beta_0.350-1.620_128_PERIODIC.tsv');


T = data(:,1);
E = data(:,2); stdE = data(:,3)/sqrt(2e7);
rhox = data(:,4); rhoy = data(:,5)/sqrt(2e7);




figure(1)
subplot(1,2,1)
%errorbar(T,E,stdE, 'b.')
plot(T,stdE.^2./T.^2,'.')
%plot(T,E, 'b.')

set(gca, 'xScale', 'lin');
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$\langle E\rangle/J$', 'interpreter', 'LaTeX');



subplot(1,2,2)
%errorbar(T,abs(M),stdM, 'b.')
%plot(T,abs(M), 'b.')
plot(T,rhox, '.'), hold on
plot(T,T);

set(gca, 'yScale', 'lin', 'ylim', [0,1])
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$\langle M\rangle$', 'interpreter', 'LaTeX');















