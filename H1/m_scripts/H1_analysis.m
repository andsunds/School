tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
set(0,'DefaultFigureWindowStyle','docked');
GRAY = 0.7*[0.9 0.9 1];

%% task 1 : lattice energies
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

ax = gca;
ax.YLim = [-13.45 -13.42];
h1 = plot( v_min*[1 1], ax.YLim, '--k'); % plot vertical line at v_min


ax.YTick = (-13.45:0.01:-13.42);
ylabel('$E_{\rm pot}$ [eV/unit cell]');
xlabel('$a_0^3$ [\AA$^3$]');
legend('data', 'quadratic fit', ['$V_{\rm eq} \approx \, $' num2str(round(v_min,2)) '\, \AA$^3$'], ...
    'location', 'southeast')
ImproveFigureCompPhys(gcf); h1.LineWidth = 2; setFigureSize(gcf);

%axis([63 68 ylim(1) 0]);
saveas(gcf, '../figures/potential_energy.eps', 'epsc')


%% task 2
%clc;
clf;clear

dt=[1e-2,5e-3,2e-3,1e-3];
for i=1:4
T_data = load(sprintf('../data/temperature_dt-%0.0e_Task2.tsv',dt(i)));
E_data = load(sprintf('../data/total_energy_dt-%0.0e_Task2.tsv',dt(i)));
t = T_data(:,1);
T = T_data(:,2);
E = E_data(:,2);

t_eq=0.5;


fprintf('dt = %0.0e\n',dt(i));

T_avg=mean(T(t>t_eq));
T_std=std(T(t>t_eq));
fprintf('\tT = %0.2f +- %0.1f %%\n', T_avg, abs(T_std/T_avg)*100);

E_avg=mean(E(t>t_eq));
E_std=std(E(t>t_eq));
fprintf('\tE = %0.2f +- %0.1e %%\n', E_avg, abs(E_std/E_avg)*100);

figure(i);clf
plot(t, T)
yyaxis right
plot(t, E)
ylim(E_avg*(1+0.001*[1,-1]));
end


%% test production pressure and temp
clc; clf;
%clear



%data = load(sprintf('../data/temperature_dt-1e-02_Task3.tsv'));
%data = load('../data/temp-700_pres-1_Task3.tsv');
data = load('../data/temp-500_pres-1_Prod-test.tsv');
%data = load('../data/temp-700_pres-1_Prod-test.tsv');
bar = 6.2415e-07;

t = data(:,1);
T = data(:,2)-273.15;
P = data(:,3)/bar;


t_eq=0.5;


%fprintf('dt = %0.0e\n',dt(i));

T_avg=mean(T(t>t_eq));
T_std=std(T(t>t_eq));
fprintf('\tT = %0.2f +- %0.1f %%\n', T_avg, abs(T_std/T_avg)*100);

P_avg=mean(P(t>t_eq));
P_std=std(P(t>t_eq));
fprintf('\tP = %0.2f +- %0.1f %%\n', P_avg, abs(P_std/P_avg)*100);

yyaxis left
plot(t,T, 'color', GRAY),hold on
plot(t, cumsum(T)./(1:length(t))','-k')
ylabel('$T \, [^\circ \rm C]$')

ylim([400,800])


yyaxis right
plot(t,P),hold on
plot(t,cumsum(P)./(1:length(t))','-k')


ylabel('$P \,[\rm bar]$')
ylim([-50,200])

xlabel('$t$\, [ps]')

ImproveFigureCompPhys(gcf, 'linewidth', 3, 'LineColor', {'MYORANGE', GRAY, 'MYBLUE', GRAY}');

%% determine displacements and MSD
temperatures = num2str([500;700]);
clc; clf;
figure(10); clf;
FILENAMES = strcat({'../data/temp-'}, temperatures, '_pres-1_displacements.tsv');
FILENAMES_Dyn = strcat({'../data/temp-'}, temperatures, '_pres-1_dynamicProperties.tsv');
FILENAMES_Pow = strcat({'../data/temp-'}, temperatures, '_pres-1_power-spectrum.tsv');
for iFile = 1:numel(FILENAMES)
    
    figure(iFile); clf;
    data = load(FILENAMES{iFile});
    t = data(:,1);
    dx = data(:,2:end);
    
    plot(t, dx.^2); hold on;
    

    data = load(FILENAMES_Dyn{iFile});
    MSD = data(:,2);
    vel_corr = data(:,3);
    plot(t, MSD, 'k')
    
    if iFile ==2 % liquid
        tStart = 1;
        D = MSD(t>tStart)./(6*t(t>tStart));
        selfDiffusionCoeff = mean(D); % in Å^2 /ps
        plot(t, 6*t*selfDiffusionCoeff, ':');
    end
    
    leg = legend( strcat({'$n=$'}, num2str((1:size(dx,2))'))');
    leg.Location='northwest';
    xlabel('$t$ [ps]')
    ylabel('$\Delta x^2 \,[\rm \AA^2]$')
    if iFile ==1
        ylim([ 0 1.0]);
    else
        ylim([0 200]);
    end
    ImproveFigureCompPhys(gcf);
    
    figure(10)
    plot(t, vel_corr/vel_corr(1)); hold on;
    xlim([0 1])
    
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
    N_times = round(length(t)/2); % we have too bad statistics at later times. 
    deltaf = 1/(N_times * dt);
    omegavec = 0:deltaf:(1/(2*dt));
    %PhiHat = 2 * trapz(t(1:N_times), (vel_corr(1:N_times) * ones(size(omegavec))) .* cos(t(1:N_times) * omegavec ), 1); %dimension 1
    PhiHat =  1/2 * 1/N_times * 2 * sum( (vel_corr(1:N_times) * ones(size(omegavec))) .* cos(t(1:N_times) * omegavec ), 1); %dimension 1
    
    figure(11); 
    
    plot(omegavec/(2*pi), PhiHat); hold on;
    plot(freq, pow_spec, ':'); hold on;
    if iFile ==2 % liquid
        tStart = 1;
        selfDiffusionCoeff_spectral = PhiHat(1)/6; % in Å^2 /ps
    end
    
end

disp([selfDiffusionCoeff selfDiffusionCoeff_spectral]);

figure(10)
xlim([0 1])
leg = legend(strcat({'$T='}, num2str([500;700]), '\,^\circ $C'));
leg.Location='northeast';
xlabel('$t$ [ps]')
ylabel('$\Phi (t)/\Phi(0)$')


figure(11)
leg = legend('$T= 500 \, ^\circ $C, $ \hat \Phi$' , '$T= 500 \, ^\circ $C, $|\hat v|^2$',...
    '$T= 700 \, ^\circ $C, $ \hat \Phi$', '$T= 700 \, ^\circ $C, $|\hat v|^2$');
xlim([0 30])
ylim([0 Inf])
ImproveFigureCompPhys('LineColor', {'r', 'MYRED', 'GERIBLUE','MYLIGHTBLUE'}');
%%
clc;clf;


FILENAME = '../data/INIDATA_temp-700_pres-1.bin';
fID=fopen(FILENAME,'rb');
data1=fread(fID,[3,inf],'real*8').';
fclose(fID);

AMU = 1.0364e-4;
m_Al = 27*AMU;
kB= 8.61733e-5;
N_atoms=4^4;

T=sum(sum(data1.^2,2),1) / (3*m_Al*N_atoms*kB)


data2=load('../data/phase-space_temp-500_pres-1.tsv');

T=sum(sum(data2(:,4:end).^2,2),1) / (3*m_Al*N_atoms*kB)




















