%% initial

tmp = matlab.desktop.editor.getActive; %% cd to current path
cd(fileparts(tmp.Filename));
set(0,'DefaultFigureWindowStyle','docked');
warning('off','MATLAB:handle_graphics:exceptions:SceneNode'); % interpreter
GRAY = 0.7*[0.9 0.9 1];
AMU = 1.0364e-4;
m_Al = 27*AMU;
%% task 1: lattice energies
clc

energy_data = load('../data/lattice_energies.tsv');
a0 = energy_data(:,1);
v0 = a0.^3;

energy = energy_data(:,2);
figure(1);clf;
plot(v0,energy, 'xk');

start_v = 64;
end_v = 68;
indToInclude = (v0 > start_v) & (v0 < end_v);
p = polyfit(v0(indToInclude),energy(indToInclude),2);
hold on;

vvec = linspace(start_v, end_v);
plot(vvec, p(1)*vvec.^2 + p(2)*vvec + p(3), '-r');
xlim([64 68]);

v_min = -p(2)/(2*p(1));
a_min = v_min^(1/3);
omega_res = sqrt(2*p(1)*a_min^4/m_Al);
f_res = omega_res/(2*pi); % rough estimation of resonance frequency (?)

ax = gca; ax.YLim = [-13.45 -13.42];
h1 = plot( v_min*[1 1], ax.YLim, '--k'); % plot vertical line at v_min

ax.YTick = (-13.45:0.01:-13.42);
ylabel('$E_{\rm pot}$ [eV/unit cell]');
xlabel('$a_0^3$ [\AA$^3$]');
legend('data', 'quadratic fit', ['$V_{\rm eq} \approx \, $' ...
    num2str(round(v_min,2)) '\, \AA$^3$'], ...
    'location', 'southeast')
ax = gca; ax.Children = ax.Children(3:-1:1);
ImproveFigureCompPhys(gcf); h1.LineWidth = 2; setFigureSize(gcf, 300, 600);
saveas(gcf, '../figures/potential_energy.eps', 'epsc')
%% task 2: find a suitable timestep
clc;

dt=[2e-2,1e-2,5e-3,2e-3]; 
t_eq=0.5;

figure(1);clf;figure(2);clf;

for i=1:numel(dt)
    T_data = load(sprintf('../data/temperature_dt-%0.0e_Task2.tsv',dt(i)));
    E_data =load(sprintf('../data/total_energy_dt-%0.0e_Task2.tsv',dt(i)));
    t = T_data(:,1);
    T = T_data(:,2);
    E = E_data(:,2);
    
    fprintf('dt = %0.0e\n',dt(i));
    
    T_avg=mean(T(t>t_eq));
    T_std=std(T(t>t_eq));
    fprintf('\tT = %0.2f +- %0.1f %%\n', T_avg, abs(T_std/T_avg)*100);
    
    E_avg=mean(E(t>t_eq));
    E_std=std(E(t>t_eq));
    fprintf('\tE = %0.2f +- %0.1e %%\n', E_avg, abs(E_std/E_avg)*100);
    
    figure(1)
    plot(t, T); hold on;
    
    figure(2)
    plot(t, E);hold on;
end
for ifig = 1:2
    figure(ifig);
    h = legend(strcat({'$dt = $ '}, num2str(round(dt',4)) , ' ps'));
    xlabel('$t$ [ps]');
    ax = gca;
    if ifig ==1
        ylabel('$\mathcal T$ [K]')
        ax.YLim = [400 1800];
    else
        ylabel('$E_{\rm tot}$ [eV/unit cell]');
        ax.YTick = (-13:0.1:-10);
        ax.YLim = [-12.75 -12.5];
    end
    ImproveFigureCompPhys(gcf,'Linewidth', 2);setFigureSize(gcf, 400, 400);
end
saveas(1, '../figures/dt-scan-temperature.eps', 'epsc')
saveas(2, '../figures/dt-scan-energy.eps', 'epsc')
%% task 3: temperature and pressure equilibration, 
% and task4: test production pressure and temperature

clc; figure(10);clf;
temps = [500 700 500 700];
temperatures_str = num2str([500;700]);
FILENAMES = [strcat({'../data/temp-'}, temperatures_str,...
    '_pres-1_Task3.tsv');
    strcat({'../data/temp-'}, temperatures_str, '_pres-1_Prod-test.tsv')];
bar = 6.2415e-07;
Kelvin_to_degC = -273.15;
t_eqs = [1 1 0.5 0.5]; % approximate equilibration time
N_average_points = 100;
dt = 5e-3;
tau_equilibration = 100*dt;

for iFile = 1:numel(FILENAMES)
    figure(iFile);clf;
    data = load(FILENAMES{iFile});
    
    t = data(:,1);
    T = data(:,2)+Kelvin_to_degC;
    P = data(:,3)/bar/1e4; %GPa
    
    
    t_eq=t_eqs(iFile);

    %fprintf('dt = %0.0e\n',dt(i));
    T_avg=mean(T(t>t_eq));
    T_std=std(T(t>t_eq));
    fprintf('\tT = %0.2f +- %0.1f K\n', T_avg, abs(T_std));
    
    P_avg=mean(P(t>t_eq));
    P_std=std(P(t>t_eq));
    fprintf('\tP = %0.2f +- %0.1f bar\n', P_avg, abs(P_std));
    
    yyaxis left
    %subplot(2,1,1)
    
    if iFile <=2 % equlibration run, otherwise production
        plot(t./tau_equilibration,T, 'color', GRAY),hold on; 
        plot(t./tau_equilibration, movmean(T,N_average_points),'-k')
    else
        plot(t,T, 'color', GRAY),hold on; 
        plot(t, cumsum(T)./(1:length(t))','-k')
    end
    ylabel('$T \, [^\circ \rm C]$')
    
    if iFile <=2 % equlibration run, otherwise production
        ylim(temps(iFile)*(1+ 0.3*[-3*0.6,0.6]))
        yyaxis right
        %legend('$\mathcal{T}$', 'mov avg');
        %subplot(2,1,2)
        plot(t./tau_equilibration,P),hold on;
        plot(t./tau_equilibration, movmean(P,N_average_points),'-k')
        
        %legend('$\mathcal{P}$', 'mov avg');
        legend('$\mathcal{T}$', 'mov avg','$\mathcal{P}$', 'mov avg', 'location', 'west');
        xlabel('$t/\tau_{\rm eq}$')
        %xlim([0 5])
    else
        ylim(temps(iFile)+ 100*[-3,1])
        yyaxis right
        %legend('$\mathcal{T}$', 'cum avg');
        %subplot(2,1,2)
        plot(t,P),hold on;
        plot(t, cumsum(P)./(1:length(t))','-k')
        legend('$\mathcal{T}$', 'cum avg','$\mathcal{P}$', 'cum avg', 'location', 'west');
        %legend('$\mathcal{P}$', 'cum avg');
        xlabel('$t$\, [ps]')
    end
    ylabel('$P \,[\rm GPa]$')
    ylim([-0.6,0.6*3])
    ImproveFigureCompPhys(gcf, 'linewidth', 3, 'LineColor', ...
        {'MYORANGE', GRAY, 'MYBLUE', GRAY}');
    %ImproveFigureCompPhys(gcf, 'linewidth', 3, 'LineColor', ...
    %    {'k', GRAY}');
    setFigureSize(gcf, 400, 400); 
    if iFile <=2 % plot a0 too
        figure(10)
        a0 = data(:,4);
        plot(t./tau_equilibration, a0);hold on;
        
    end
end

saveas(1, '../figures/TP-eq-500.eps', 'epsc')
saveas(2, '../figures/TP-eq-700.eps', 'epsc')
saveas(3, '../figures/TP-prod-500.eps', 'epsc')
saveas(4, '../figures/TP-prod-700.eps', 'epsc')
saveas(10, '../figures/a0.eps', 'epsc')
figure(10);
xlabel('$t/\tau_{\rm eq}$');
ylabel('$a_0$');
legend('$T = 500^\circ$C','$T = 700^\circ$C')
setFigureSize(gcf, 300, 600);
ImproveFigureCompPhys(gcf);

%% determine displacements and MSD
temperatures_str = num2str([500;700]);
clc; clf;
figure(10); clf;
FILENAMES = strcat({'../data/temp-'}, temperatures_str, ...
    '_pres-1_displacements.tsv');
FILENAMES_Dyn = strcat({'../data/temp-'}, temperatures_str, ...
    '_pres-1_dynamicProperties.tsv');
FILENAMES_Pow = strcat({'../data/temp-'}, temperatures_str, ...
    '_pres-1_power-spectrum.tsv');
for iFile = 1:numel(FILENAMES)
    figure(iFile); clf;
    data = load(FILENAMES{iFile});
    t = data(:,1);
    dx = data(:,2:end);

    data = load(FILENAMES_Dyn{iFile});
    MSD = data(:,2);
    vel_corr = data(:,3);
    plot(t, MSD, 'k'); hold on;
    
    if iFile ==2 % liquid
        tStart = 1; tEnd = t(end)-tStart;
        goodInd = (t>tStart & t<tEnd);
        
        D = MSD(goodInd)./(6*t(goodInd));
        
        selfDiffusionCoeff = mean(D); % in Å^2 /ps
        plot(t, 6*t*selfDiffusionCoeff, ':r');
    end
    
    plot(t, dx.^2, 'color', GRAY); hold on;
    
    xlabel('$t$ [ps]')
    ylabel('$\Delta x^2 \,[\rm \AA^2]$')
    if iFile ==1
        ylim([ 0 1.0]);
        leg = legend( '$\Delta_{\rm MSD}$', 'individual trajectories');
    else
        ylim([0 50]);
        leg = legend('$\Delta_{\rm MSD}$', '$6 t D_s$', ...
            'individual trajectories');
    end
    %xlim([0 10])
 
    leg.Location='northwest';
    ImproveFigureCompPhys(gcf, 'Linewidth', 2);
    ax = gca; [ax.Children(6:end).LineWidth] = deal(5);
    ax.Children = ax.Children([6:end 1:5]);
    setFigureSize(gcf, 400, 400);
end

% velocity correlation
figure(10);clf; figure(11);clf;
n_average_points = 1;%30;
for iFile = 1:numel(FILENAMES)
    data = load(FILENAMES_Dyn{iFile});
    t = data(:,1);
    vel_corr = data(:,3);
    
    data = load(FILENAMES_Pow{iFile});
    freq = data(:,1);
    pow_spec = data(:,2);
    
    figure(10);
    plot(t, vel_corr/vel_corr(1)); hold on;
    
    dt = t(2)-t(1);
    N_times = round(length(t)/1.5);% we have bad statistics at later times.
    deltaf = 1/(N_times * dt);
    freqvec = 0:deltaf:(1/(2*dt));
    PhiHat = 2 * trapz(t(1:N_times), ...
        (vel_corr(1:N_times) * ones(size(freqvec))) .* ...
        cos(2*pi*t(1:N_times) * freqvec ), 1); 
    
    figure(11);
    plot(freqvec, PhiHat); hold on;
    plot(freq, pow_spec*t(end), '-.'); hold on;
    if iFile ==2 % liquid
        tStart = 1;
        selfDiffusionCoeff_spectral = PhiHat(1)/6; % in Å^2 /ps
    end
    
end

disp([selfDiffusionCoeff selfDiffusionCoeff_spectral]);

figure(10)
xlim([0 1]);
leg = legend(strcat({'$T='}, num2str([500;700]), '\,^\circ $C'));
leg.Location='northeast';
xlabel('$t$ [ps]')
ylabel('$\Phi (t)/\Phi(0)$')
ImproveFigureCompPhys(gcf,'LineColor', {'MYRED', 'MYLIGHTBLUE'}');
setFigureSize(gcf, 400, 400);

figure(11)
leg = legend('$T= 500 \, ^\circ $C, $ \hat \Phi$' ,...
    '$T= 500 \, ^\circ $C, $\hat P$',...
    '$T= 700 \, ^\circ $C, $ \hat \Phi$', ...
    '$T= 700 \, ^\circ $C, $\hat P$');
xlim([0 20])
ylim([0 Inf])
xlabel('$f$ [THz]')
ylabel('$\hat P$ [\AA$^2$/ps]')
setFigureSize(gcf, 300, 600);

ImproveFigureCompPhys(gcf,'LineColor', ...
    {'r', 'MYRED', 'GERIBLUE','MYLIGHTBLUE'}', 'linewidth', 3);
saveas(1, '../figures/MSD-500.eps', 'epsc')
saveas(2, '../figures/MSD-700.eps', 'epsc')
saveas(10, '../figures/Phi-t.eps', 'epsc')
saveas(11, '../figures/P-freq.eps', 'epsc')