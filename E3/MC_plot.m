data = load('MC_analysis.tsv');

k = data(:,1);
block_size = k+1;
phi = data(:,2);
VarF_norm = data(:,3);

subplot(2,1,1)
plot(k, log(phi));
ylim([-6 0]);

subplot(2,1,2);
plot(block_size, VarF_norm); hold on;
plot(block_size, VarF_norm); 
