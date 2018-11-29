%% initial

tmp = matlab.desktop.editor.getActive; %% cd to current path
cd(fileparts(tmp.Filename));
set(0,'DefaultFigureWindowStyle','docked');
warning('off','MATLAB:handle_graphics:exceptions:SceneNode'); % interpreter
GRAY = 0.7*[0.9 0.9 1];

%% task 1: MFT
figure(1);clf;
Pmin = 0;
Pmax = 1;
F_func = @(P,Tbar) -2 * P.^2 + Tbar*((1+P).*log(1+P)+(1-P).*log(1-P));
dFdP = @(P, Tbar)( -4*P + Tbar * (log(1+P)-log(1-P)) );
Tbar = linspace(0,3,100);
Peq = zeros(size(Tbar));
Pvec = linspace(-1,1,1000);
for iT = 1:numel(Tbar)
%    plot(Pvec, dFdP(Pvec, Tbar(iT))); hold on; ylim([-5 5])
    Peq(iT) = fminbnd(@(P)F_func(P, Tbar(iT)), Pmin, Pmax);
end
plot(Tbar, Peq);

xlabel('$k_B T/ \Delta E$')
ylabel('$P_{\rm equilibrium}$')
ImproveFigureCompPhys(gcf)

%% task 2: ...
clc;

Ts=[400:10:590]; 
t_eq=0;

figure(1);clf;figure(2);clf;

for i=1:numel(Ts);
    data = load(sprintf('../data/E_equilibration-T%d.tsv',Ts(i)));
    E = data(:,1);
    
    %fprintf('T = %0.0e\n',Ts(i));
    
    %E_avg=mean(E(t>t_eq));
    %E_std=std(E(t>t_eq));
    %fprintf('\tE = %0.2f +- %0.1e %%\n', E_avg, abs(E_std/E_avg)*100);
    
    
    plot(E); hold on;
end

data = load('../data/E_mean.tsv');
T = data(:,1);
Emean = data(:,2);
E_sq_mean = data(:,3);

figure(2);clf;
plot(T, Emean);
ImproveFigureCompPhys()
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
