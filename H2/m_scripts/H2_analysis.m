%% initial

tmp = matlab.desktop.editor.getActive; %% cd to current path
cd(fileparts(tmp.Filename));
set(0,'DefaultFigureWindowStyle','docked');
warning('off','MATLAB:handle_graphics:exceptions:SceneNode'); % interpreter
GRAY = 0.7*[0.9 0.9 1];
kB = 8.61733e-5;
%% task 1: MFT
doSave = 0;
clc

Pmin = 0;
Pmax = 1;

E_CuCU = -.436;
E_ZnZn = -.133;
E_CuZn = -.294;

E0=2*(E_CuCU+E_ZnZn+2*E_CuZn);
Delta_E=(E_CuCU+E_ZnZn-2*E_CuZn);

E0_bar=E0/Delta_E;
E_MFT=@(P) E0 - 2*P.^2*Delta_E;
E_MFT_bar=@(P) E0_bar - 2*P.^2;
dE_MFTdP =@(P) - 4*P*Delta_E;

F_MFT = @(P,Tbar) E_MFT_bar(P) + Tbar*(-2*log(2) + (1+P).*log(1+P)+(1-P).*log(1-P));
P_eq=@(Tbar)  fminbnd(@(P)F_MFT(P, Tbar), Pmin, Pmax, optimset('TolX',1e-9));

Tbar = linspace(0,3,1000)';
T_MFT=Tbar*Delta_E/kB;
T_MFT_degC = T_MFT - 273.15;
Peq = zeros(size(Tbar));
for iT = 1:numel(Tbar)
    Peq(iT) = P_eq(Tbar(iT));
end

% plot P(T) and make a fit
figure(1);clf
plot(Tbar, Peq);hold on

dT=2-Tbar(Tbar<2);
Peq_nonzero = Peq(Tbar<2);

I_good = (dT<0.1);
log_dT = log(dT(I_good));
log_P  = log(Peq_nonzero(I_good));
A=[ones(size(log_dT)), log_dT]\log_P;
b     = exp(A(1));
alpha = A(2);
fprintf('alpha = %.3f\n', alpha)

P_approx = @(alpha,b,Tbar) b*(2-Tbar).^alpha;
plot(Tbar(Tbar<2),P_approx(alpha,b,Tbar(Tbar<2)),'k:')
xlabel('$k_B T/ \Delta E$')
ylabel('$P$')
legend('$P$', 'fit $P \propto (2-\bar T)^\alpha$')
ylim([0 1.3]);
if doSave; setFigureSize(gcf, 300, 600); end

% plot E_MFT and the fit
figure(2);clf
plot(Tbar,E_MFT(Peq)); hold on
plot(Tbar,E_MFT(P_approx(alpha,b,Tbar)),'k:')
xlabel('$k_B T/ \Delta E$')
ylabel('$U$ [eV/cell]')
legend('$E_{\rm MFT}$', 'fit $P \propto (2-\bar T)^\alpha$', 'location', 'NorthWest');
ylim([-2.36 -2.3]);
if doSave; setFigureSize(gcf, 300, 600); end

figure(3);clf
C=diff(E_MFT(Peq))./diff(T_MFT);
plot(Tbar(1:end-1), C*1e3); hold on
C_approx=4*b^2*kB*alpha*(2-Tbar).^(2*alpha-1);
plot(Tbar(Tbar<2),1e3*C_approx(Tbar<2),'k:')
xlabel('$k_B T/ \Delta E$')
ylabel('$C$ [meV K$^{-1}$/cell]')
legend('$C_{\rm MFT}$', 'fit $P \propto (2-\bar T)^\alpha$', 'location', 'NorthWest');
ylim([0 0.3])
if doSave; setFigureSize(gcf, 300, 600); end

ImproveFigureCompPhys()
if doSave
    saveas(1, '../figures/P_MFT.eps', 'epsc');
    saveas(2, '../figures/E_MFT.eps', 'epsc');
    saveas(3, '../figures/C_MFT.eps', 'epsc');
end
%% task 2: ...
clc;
doSave = 1;
Ts=[-200:20:600]';
TsToPlot = [300 440 600]';
t_eq=0;

figure(10);clf;

for i=1:numel(TsToPlot)
    data = load(sprintf('../data/E_equilibration-T%d.tsv',TsToPlot(i)));
    E = data(:,1);
    steps = 1:length(E);
    %P = data(:,2);
    plot(steps/1e6, E*1000); hold on
end
legstr = strcat({'$T='}, num2str(TsToPlot), '^\circ$ C');
legend(legstr, 'location', 'NorthWest');
ylabel('$E$ [meV/$N_{\rm bonds}$]')
xlabel('$N_{\rm steps}/10^6$')
if doSave
    ImproveFigureCompPhys(gcf)
    setFigureSize(gcf, 300, 600); 
    saveas(gcf, '../figures/equilibration.eps', 'epsc');
end


figure(1000); clf;
[ns_Phi,ns_block] = deal(nan(size(Ts)));
Nskip = 10;
for i=1:numel(Ts)
    data = load(sprintf('../data/stat_inefficiency-T%d.tsv',Ts(i)));
    k = data(:,1);
    block_size = k+Nskip;
    phi = data(:,2);
    VarF_norm = data(:,3);
    kstar = k(find(log(phi)<-2, 1, 'first'));
    if ~isempty(kstar)
        ns_Phi(i) = kstar;
    end
    N_avg = 20;
    filtereddata = movmean(VarF_norm,N_avg);
    ns_block(i) = filtereddata(end);
    
    if any(Ts(i) == TsToPlot)
        subplot(2,1,1)
        plot(k, log(phi));hold on;
        
        plot([0 kstar kstar], [-2 -2 -6],':k')
        ylim([-4 0]);
        legend('data', 'estimated $n_s$', 'location', 'northeast');
        xlabel('$k$'); ylabel('ln $\phi_k$');
        xlim([0 2e5])
        
        subplot(2,1,2);
        plot(block_size, VarF_norm); hold on;
        
        plot(block_size(N_avg:end), filtereddata(N_avg:end));
        plot(block_size, filtereddata(end)*ones(size(block_size)), ':k');
        legend('data', 'moving average', 'estimated $n_s$', 'location', 'northwest');
        xlabel('block size $B$'); ylabel('$B$ Var[$F$]/Var[$f$] ');
        ylim([0 2e5])
    end
end
%Ts = Ts(~isnan(ns_Phi));
%ns_Phi = ns_Phi(~isnan(ns_Phi));
%ns_block = ns_block(~isnan(ns_Phi));

ImproveFigureCompPhys()
%%

data = load('../data/E_production.tsv');
T_degC = data(:,1);
N_Cu = 1e3;
N_timeSteps = 1e7;

Emean_approx = data(:,2);
Emean_shifted = data(:,3);
E_sq_mean_shifted = data(:,4);

E_Var = (E_sq_mean_shifted - Emean_shifted.^2);

Cv = 1./(kB * (T_degC+273.15).^2).*E_Var;
U = Emean_shifted + Emean_approx;
U_std = sqrt(E_Var/N_timeSteps);
P = data(:,5);
P_std = sqrt((data(:,6)-P.^2)/N_timeSteps); % without ns so far
r = data(:,7);
r_std = sqrt((data(:,8)- r.^2)/N_timeSteps);

ind = zeros(size(Ts));
for i = 1:numel(Ts)
    ind(i) = find(Ts(i) == T_degC);
end

figure(11);clf;

errorbar(Ts, U(ind), 2*U_std(ind).*sqrt(ns_Phi), '.k','linewidth', 1.5); hold on;
plot(T_degC, U); hold on;

plot(T_degC, cumtrapz(T_degC, Cv) + U(1));

figure(12); clf;
plot(T_degC, Cv/N_Cu); hold on;
plot(T_MFT_degC(1:end-1), C); hold on

figure(13);clf;
errorbar(Ts, P(ind), 2*P_std(ind).*sqrt(ns_Phi), '.k', 'linewidth', 1.5); hold on;
%errorbar(Ts, P(ind), 2*P_std(ind).*sqrt(ns_block), '.r','linewidth', 1.5);hold on;
plot(T_degC, P, 'color', GRAY); hold on;

plot(T_MFT_degC, Peq, '--k');

figure(14);clf;
errorbar(Ts, r(ind), 2*r_std(ind).*sqrt(ns_Phi), '.k','linewidth', 1.5);hold on;
hold on; plot(T_degC, r, T_degC, P.^2, T_MFT_degC, Peq.^2, 'k');

legend('$r$','$P^2$',  '$r_{\rm MFT}$ ')
ImproveFigureCompPhys('linewidth', 2)

% for ifig = 1:2
%     figure(ifig);
%     h = legend(strcat({'$dt = $ '}, num2str(round(dt',4)) , ' ps'));
%     xlabel('$t$ [ps]');
%     ax = gca;
%     if ifig ==1
%         ylabel('$T$ [K]')
%         ax.YLim = [400 1800];
%     else
%         ylabel('$E_{\rm tot}$ [eV/unit cell]');
%         ax.YTick = (-13:0.1:-10);
%         ax.YLim = [-12.6 -12.0];
%     end
%     ImproveFigureCompPhys(gcf,'Linewidth', 2);setFigureSize(gcf, 400, 400);
% end
% saveas(1, '../figures/dt-scan-temperature.eps', 'epsc')
% saveas(2, '../figures/dt-scan-energy.eps', 'epsc')

%%


