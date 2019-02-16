addpath('../H1/m_scripts/');

data = load('MC_analysis.tsv');

k = data(:,1);
block_size = k+1;
phi = data(:,2);
VarF_norm = data(:,3);

clf;
subplot(2,1,1)
plot(k, log(phi));
kstar = find(log(phi)<-2, 1, 'first');hold on;
plot([0 kstar kstar], [-2 -2 -6],':k') 
ylim([-6 0]);
legend('data', 'estimated $n_s$', 'location', 'northeast');
xlabel('$k$'); ylabel('$\phi_k$'); 

subplot(2,1,2);
plot(block_size, VarF_norm); hold on;
filtereddata = movmean(VarF_norm,50);
plot(block_size, filtereddata); 
plot(block_size, filtereddata(end)*ones(size(block_size)), ':k'); 
legend('data', 'moving average', 'estimated $n_s$', 'location', 'southeast');
xlabel('block size $B$'); ylabel('$B$ Var[$F$]/Var[$f$] '); 

ImproveFigureCompPhys();

saveas(gcf, 'MC_plot.pdf')