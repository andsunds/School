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
legend('$U_{\rm MFT}$', 'fit $P \propto (2-\bar T)^\alpha$', 'location', 'NorthWest');
ylim([-2.36 -2.3]);
if doSave; setFigureSize(gcf, 300, 600); end

figure(3);clf
C_MFT=diff(E_MFT(Peq))./diff(T_MFT);
plot(Tbar(1:end-1), C_MFT*1e3); hold on
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


%% task 2: equilibration and statistical inefficiency
clc;
doSave = 0;
Ts=[-200:20:600]';
TsToPlot = [300 440 600]';
t_eq=0;

figure(1);clf;

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
ImproveFigureCompPhys(1)

figure(3); clf;figure(2); clf;
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
    N_avg = 100;
    filtereddata = movmean(VarF_norm,N_avg);
    ns_block(i) = filtereddata(end);
    
    if any(Ts(i) == TsToPlot)
        figure(2)
        
        semilogx(k, log(phi));hold on;
        plot([0.1 kstar kstar], [-2 -2 -6],':k')
        
        figure(3)
        semilogy(block_size, VarF_norm, '.'); hold on;
        plot(block_size(N_avg:end), filtereddata(N_avg:end));
        plot(block_size, filtereddata(end)*ones(size(block_size)), ':k');
        
    end
end

figure(4); clf;
plot(Ts, ns_Phi, 'k',Ts, ns_block, '--r')
ax = gca; 
ax.YTickLabel = {'0', '$10^5$', '$2\cdot 10^5$','$3\cdot 10^5$','$4\cdot 10^5$','$5\cdot 10^5$'}';
ylabel('$n_s$');
legend('correlation function $\Phi$', 'block average');
xlabel('$T$ [$^\circ$C]');
ImproveFigureCompPhys(gcf)

legs_Phi = cell(6,1);
legs_block = cell(9,1);
for i = 1:numel(TsToPlot)
    tt = ['$T=' num2str(TsToPlot(i)) '$ K: '];
    legs_Phi{1 + 2*(i-1)} = [tt 'data'];
    legs_Phi{2 + 2*(i-1)} = 'estimated $n_s$';
    legs_block{1 + 3*(i-1)} = [tt 'data'];
    legs_block{2 + 3*(i-1)} = 'moving average';
    legs_block{3 + 3*(i-1)} = 'estimated $n_s$';
end

figure(2);

legend(legs_Phi, 'location', 'northeastoutside');
xlabel('$k$'); ylabel('ln $\phi_k$');
ylim([-3.5 0]);
xlim([2e3 3e5])
%ax = gca; ax.XTick = [3e3 1e4 3e4 1e5 3e5];
%ax.XTickLabel = {'$3\cdot 10^3$', '$10^4$','$3\cdot 10^4$','$10^5$','$3\cdot 10^5$'}';
figure(3);
ax = gca;
[ax.Children(:).MarkerSize] = deal(12);
legend(legs_block, 'location', 'northeastOutSide');
xlabel('block size $B$'); 
ylabel('$B$ Var[$F$]/Var[$f$]');
ylim([2e3 2e5])
ax = gca; 
ax.XTickLabel = {'0', '$10^5$', '$2\cdot 10^5$','$3\cdot 10^5$','$4\cdot 10^5$','$5\cdot 10^5$'}';

ImproveFigureCompPhys(2, 'LineColor', {'LINNEAGREEN','LINNEAGREEN','GERIBLUE','GERIBLUE', 'k', 'k'}',...
    'LineStyle', {':','-.',':','-', ':', '--'}')
ImproveFigureCompPhys(3, 'LineColor', {'LINNEAGREEN','LINNEAGREEN','LINNEAGREEN',...
    'GERIBLUE','GERIBLUE','GERIBLUE', 'k', 'k', 'k'}',...
    'LineStyle', {':','-.','none',':','-','none', ':', '--','none'}');
if doSave
    figure(1);
    setFigureSize(gcf, 300, 600); 
    saveas(gcf, '../figures/equilibration.eps', 'epsc');
    figure(2);
    setFigureSize(gcf, 350, 900); 
    saveas(gcf, '../figures/stat_inefficiency_Phi.eps', 'epsc');
    figure(3);
    setFigureSize(gcf, 350, 900); 
    saveas(gcf, '../figures/stat_inefficiency_block.eps', 'epsc');
    figure(4);
    setFigureSize(gcf, 300, 600); 
    saveas(gcf, '../figures/stat_inefficiency_both.eps', 'epsc');
end


%% task 2: U, C, P and r

doSave = 0;

data = load('../data/E_production.tsv');
T_degC = data(:,1);
N_Cu = 1e3;
N_timeSteps = 1e7;

Emean_approx = data(:,2)/N_Cu; % divide by N_Cu to get energy and Cv per cell
Emean_shifted = data(:,3)/N_Cu;
E_sq_mean_shifted = data(:,4)/N_Cu^2; 

E_Var = (E_sq_mean_shifted - Emean_shifted.^2);

Cv = 1./(kB * (T_degC+273.15).^2).*E_Var*N_Cu;
U = (Emean_shifted + Emean_approx); 
U_std = sqrt(E_Var/N_timeSteps);
P = data(:,5);
P_std = sqrt((data(:,6)-P.^2)/N_timeSteps); % without ns so far
r = data(:,7);
r_std = sqrt((data(:,8)- r.^2)/N_timeSteps);

ind = zeros(size(Ts));
for i = 1:numel(Ts)
    ind(i) = find(Ts(i) == T_degC);
end

figure(1);clf;
plot(T_degC, U); hold on;
errorbar(Ts, U(ind), 2*U_std(ind).*sqrt(ns_Phi), '.k','linewidth', 2.5); hold on;
plot(T_MFT_degC, E_MFT(Peq), '-.'); hold on
ImproveFigureCompPhys(gcf, 'LineColor', {'GERIBLUE', 'r'}');
legend('$U$', '$U\pm 2 \sigma$ (with $n_{s, \rm \Phi})$', '$E_{\rm MFT}$', 'Location', 'NorthWest');
ylabel('$U$ [eV/cell]')

figure(2); clf;
plot(T_degC(2:end), 1e3*diff(U)./diff(T_degC)); hold on;
plot(T_degC, 1e3*Cv); 
plot(T_MFT_degC(1:end-1), 1e3*C_MFT, '-.');
ImproveFigureCompPhys(gcf, 'LineColor', {'GERIBLUE', 'k',GRAY}');
legend('$C, {\rm Var(E)}$', '$C, {\partial U/ \partial T}$', '$C_{\rm MFT}$', 'Location', 'NorthWest');
ylabel('$C$ [meV/cell]')

figure(3);clf;
plot(T_degC, P, 'r'); hold on;
errorbar(Ts, P(ind), 2*P_std(ind).*sqrt(ns_Phi), '.k', 'linewidth', 2.5); hold on;
plot(T_MFT_degC, Peq, '-.k');
ImproveFigureCompPhys(gcf, 'LineColor', {'GERIBLUE', 'r'}');
legend('$P$', '$P\pm 2 \sigma$ (with $n_{s, \rm \Phi})$', '$P_{\rm MFT}$', 'Location', 'SouthWest');
ylabel('$P$ ')


figure(4);clf;
plot(T_degC, r, 'r');hold on;
errorbar(Ts, r(ind), 2*r_std(ind).*sqrt(ns_Phi), '.k','linewidth', 1.5);hold on;
plot(T_degC, P.^2, '--',T_MFT_degC, Peq.^2, '-.');
ImproveFigureCompPhys(gcf, 'LineColor', {'GERIBLUE', 'LINNEAGREEN','r'}');
legend('$r$', '$r\pm 2 \sigma$ (with $n_{s, \rm \Phi})$', '$P^2$','$r_{\rm MFT}$', 'Location', 'SouthWest');
ylabel('$r$ ')

ImproveFigureCompPhys((2:4), 'linewidth', 2)

if doSave
    for ifig = 1:4;
        figure(ifig)
        setFigureSize(gcf, 300, 600); 
        xlabel('$T$ [$^\circ$C]');
        axis tight
        xlim([-200 Inf])
    end
    ImproveFigureCompPhys(1:4);
    saveas(1, '../figures/U.eps', 'epsc');
    saveas(2, '../figures/C.eps', 'epsc');
    saveas(3, '../figures/P.eps', 'epsc');
    saveas(4, '../figures/r.eps', 'epsc');
end

