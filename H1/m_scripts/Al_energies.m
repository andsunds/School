tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath('../../../../scripts');
set(0,'DefaultFigureWindowStyle','docked');
clc

energy_data = load('../data/lattice_energies.tsv');
a0 = energy_data(:,1);
v0 = a0.^3;

energy = energy_data(:,2);
figure(1);clf;
plot(v0,energy, 'x-', 'color', 0.7*[0.9 0.9 1]);
%ImproveFigure(gcf, 'profile', 'paper1'); ImproveFigure(gcf, 'profile', 'paper1')
ylabel('$E_{\rm pot}$ [eV/unit cell]');
xlabel('$a_0^3$ [\AA$^3$]');


start_a = 4;
end_a = 68.^(1/3);
indToInclude = (a0 > start_a) & (a0 < end_a);
p = polyfit(v0(indToInclude),energy(indToInclude),2);
hold on;

vvec = linspace(start_a^3, end_a^3);
plot(vvec, p(1)*vvec.^2 + p(2)*vvec + p(3), '-r');
xlim([64 68]);

v_min = -p(2)/(2*p(1));
ax = gca;
%ylim = ax.YLim;
ylim([-13.45 -13.42]);
h1 = plot( v_min*[1 1], ylim, '--k');
legend('data', 'quadratic fit', ['$V_{\rm min} \approx \, $' num2str(round(v_min,3))], ...
    'location', 'southeast')

%axis([63 68 ylim(1) 0]);
%ImproveFigure(gcf, 'profile', 'paper1')
%h1.LineWidth = 2;


%%
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


%%
clc; clf;clear


%data = load(sprintf('../data/temperature_dt-1e-02_Task3.tsv'));
%data = load('../data/temp-500_pres-1_Task3.tsv');
data = load('../data/temp-500_pres-1_Prod-test.tsv');
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

%figure(i);clf
plot(t,T),hold on
%plot(t, cumsum(T)./(1:length(t))','-k')
%ylim([0,1000])
yyaxis right
plot(t,P),hold on
%plot(t,cumsum(P)./(1:length(t))','-k')
%plot(t,1+500*exp(-t/0.5))
%ylim([-100,200])

%%
clc;clf;clear


FILENAME = '../data/mom_temp-500_pres-1.bin';
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




















