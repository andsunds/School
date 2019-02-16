

%%
tmp = matlab.desktop.editor.getActive; %% cd to current path
cd(fileparts(tmp.Filename));
set(0,'DefaultFigureWindowStyle','docked');
warning('off','MATLAB:handle_graphics:exceptions:SceneNode'); % interpreter
%%

viscosity = 18.6e-6;

tau = 48.5e-6;

rho = 2.65e3;
D = 2.79e-6;
r = D/2;
m = rho * 4 * pi * r^3/3;

tau_Stokes = m/(6*pi*viscosity * r)

% 61 vs 48 = consistent. 

%% tau small

data = load('xv_tau-48.tsv');
t = data(:,1);
x = data(:, 2:2:10);
v = data(:, 3:2:11);
data = load('mu_sigma_tau-48.tsv');
mux = data(:,2);
sigmax = data(:,3);
muv = data(:,4);
sigmav = data(:,5);

figure(1);clf;
plot(t, 1e3*x,'linewidth', 1.4); hold on;
plot(t, 1e3*mux, '-k', t, 1e3*(mux-sigmax), ':k', t, 1e3*(mux+sigmax), ':k', 'linewidth', 3)
xlabel('t [\mu s]')
ylabel('x [n m]')

figure(2);clf;
plot(t, 1e3*v,'linewidth', 1.4); hold on;
plot(t, 1e3*muv, '-k', t, 1e3*(muv-sigmav), ':k', t, 1e3*(muv+sigmav), ':k', 'linewidth', 3)
xlabel('t [\mu s]')
ylabel('v [\mu m/ms]')

data = load('dist_tau-48.tsv');
rho_t = data(1,:);
rho_x = data(2:end, 1:2:end);
rho_v = data(2:end, 2:2:end);
xlim = 0.15;
vlim = 2e-3;
N_dist_points = size(data,1)-1;
dx = 2*xlim / N_dist_points;
dv = 2*vlim / N_dist_points;
xvec = -xlim:dx:(xlim-dx);
vvec = -vlim:dv:(vlim-dv);
figure(3);clf;
plot(xvec', rho_x);
legend(strcat({'t='}, num2str(rho_t(1:2:end)'), '\mu s'), 'location', 'northwest')
xlabel('\mu m')
ylabel('\rho_x [#particles/\mu m]' )
ylim([0 50])

figure(4);clf;
plot(vvec'*1e3, rho_v*1e-3);
legend(strcat({'t='}, num2str(rho_t(1:2:end)'),'\mu s'),'location', 'northwest')
xlabel('\mu m/ms')
ylabel('\rho_v [#particles/ (\mu m/ms)]')
ylim([0 2])

saveas(1, 'tausmall_x.eps', 'epsc');
saveas(2, 'tausmall_v.eps', 'epsc');
saveas(3, 'tausmall_rhox.eps', 'epsc');
saveas(4, 'tausmall_rhov.eps', 'epsc');
%% tau large

data = load('xv_tau-147.tsv');
t = data(:,1);
x = data(:, 2:2:10);
v = data(:, 3:2:11);
data = load('mu_sigma_tau-147.tsv');
mux = data(:,2);
sigmax = data(:,3);
muv = data(:,4);
sigmav = data(:,5);

figure(1);clf;
plot(t, 1e3*x,'linewidth', 1.4); hold on;
plot(t, 1e3*mux, '-k', t, 1e3*(mux-sigmax), ':k', t, 1e3*(mux+sigmax), ':k', 'linewidth', 3)
xlabel('t [\mu s]')
ylabel('x [n m]')

figure(2);clf;
plot(t, 1e3*v,'linewidth', 1.4); hold on;
plot(t, 1e3*muv, '-k', t, 1e3*(muv-sigmav), ':k', t, 1e3*(muv+sigmav), ':k', 'linewidth', 3)
xlabel('t [\mu s]')
ylabel('v [\mu m/ms]')

data = load('dist_tau-147.tsv');
rho_t = data(1,:);
rho_x = data(2:end, 1:2:end);
rho_v = data(2:end, 2:2:end);
xlim = 0.15;
vlim = 2e-3;
N_dist_points = size(data,1)-1;
dx = 2*xlim / N_dist_points;
dv = 2*vlim / N_dist_points;
xvec = -xlim:dx:(xlim-dx);
vvec = -vlim:dv:(vlim-dv);
figure(3);clf;
plot(xvec', rho_x);
legend(strcat({'t='}, num2str(rho_t(1:2:end)'), '\mu s'), 'location', 'northwest')
xlabel('\mu m')
ylabel('\rho_x [#particles/\mu m]' )
ylim([0 50])

figure(4);clf;
plot(vvec'*1e3, rho_v*1e-3);
legend(strcat({'t='}, num2str(rho_t(1:2:end)'),'\mu s'),'location', 'northwest')
xlabel('\mu m/ms')
ylabel('\rho_v [#particles/ (\mu m/ms)]')
ylim([0 2])
saveas(1, 'taularge_x.eps', 'epsc');
saveas(2, 'taularge_v.eps', 'epsc');
saveas(3, 'taularge_rhox.eps', 'epsc');
saveas(4, 'taularge_rhov.eps', 'epsc');