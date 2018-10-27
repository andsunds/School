%% 4b
clc;clear;clf

t=linspace(0,2).';

c=(-2:8);
C=9/4;

x = repmat(c, length(t),1)-repmat((t+.5).^2,1,length(c));
X = C-(t+.5).^2;

plot(x,t, ':k'); hold on
plot(X,t, '-k', 'linewidth',3)
plot(0,1, 'ok', 'markersize', 10)
set(gca, 'ytick', [0:.5:2], 'fontsize',15);
grid on;
axis([-3,3,0,2]);

xlabel('$x$', 'interpreter', 'latex')
ylabel('$t$', 'interpreter', 'latex')