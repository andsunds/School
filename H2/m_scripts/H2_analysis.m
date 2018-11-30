%% initial

tmp = matlab.desktop.editor.getActive; %% cd to current path
cd(fileparts(tmp.Filename));
set(0,'DefaultFigureWindowStyle','docked');
warning('off','MATLAB:handle_graphics:exceptions:SceneNode'); % interpreter
GRAY = 0.7*[0.9 0.9 1];

%% task 1: MFT
clc

Pmin = 0;
Pmax = 1;

kB= 8.61733e-5;

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
%dFdP = @(P, Tbar)( -4*P + Tbar * (log(1+P)-log(1-P)) );

P_eq=@(Tbar)  fminbnd(@(P)F_MFT(P, Tbar), Pmin, Pmax, optimset('TolX',1e-9));

Tbar = linspace(0,3,1000)';
T=Tbar*Delta_E/kB;
Peq = zeros(size(Tbar));
for iT = 1:numel(Tbar)
    Peq(iT) = P_eq(Tbar(iT));
end

figure(1);clf
plot(T, Peq);

xlabel('$k_B T/ \Delta E$')
ylabel('$P_{\rm equilibrium}$')


figure(2);clf

dT=2-Tbar(Tbar<2);
Peq_nonzero = Peq(Tbar<2);

plot(dT, Peq_nonzero); hold on
set(gca,'xscale','log','yscale','log');

I_good = (dT<0.1);
log_dT = log(dT(I_good));
log_P  = log(Peq_nonzero(I_good));

A=[ones(size(log_dT)), log_dT]\log_P;
b     = exp(A(1));
alpha = A(2);
fprintf('alpha = %.3f\n', alpha)

P_approx = @(alpha,b,Tbar) b*(2-Tbar).^alpha;

plot(dT,P_approx(alpha,b,2-dT),'k:')



xlabel('$k_B T/ \Delta E$')
ylabel('$P_{\rm equilibrium}$')

figure(3);clf
plot(T,E_MFT(Peq)); hold on

plot(T,E_MFT(P_approx(alpha,b,Tbar)),'k:')


figure(4);clf
%dPdT=diff(Peq)./diff(T);
%dE_MFTdP(Peq(1:end-1))*
C=diff(E_MFT(Peq))./diff(T);
plot(Tbar(1:end-1), C); hold on

C_approx=4*b^2*kB*alpha*(2-Tbar).^(2*alpha-1);
plot(Tbar(Tbar<2),C_approx(Tbar<2),'k:')


ImproveFigureCompPhys()

%% task 2: ...
clc;

Ts=[300:50:600]; 
%Ts=[300:50:450]; 
t_eq=0;

figure(10);clf;

for i=1:numel(Ts)
    data = load(sprintf('../data/E_equilibration-T%d.tsv',Ts(i)));
    E = data(:,1);
    P = data(:,2);
    
    %fprintf('T = %0.0e\n',Ts(i));
    
    %E_avg=mean(E(t>t_eq));
    %E_std=std(E(t>t_eq));
    %fprintf('\tE = %0.2f +- %0.1e %%\n', E_avg, abs(E_std/E_avg)*100);
    
    %plot(P); hold on
    plot(E); hold on;
end

data = load('../data/E_mean.tsv');
T = data(:,1);
Emean = data(:,2);
E_sq_mean = data(:,3);
%%
figure(11);clf;
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
