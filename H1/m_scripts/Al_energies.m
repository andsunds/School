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
ImproveFigure(gcf, 'profile', 'paper1'); ImproveFigure(gcf, 'profile', 'paper1')
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
ImproveFigure(gcf, 'profile', 'paper1')
%h1.LineWidth = 2;


%%
tmp_data = load('../data/atom_energies_tmp.tsv');
t = tmp_data(:,1);
energy = tmp_data(:,2);

clf;
plot(t, energy)

