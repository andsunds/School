% plot the displacements
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('disp.dat');

%plot
figure;
plot(data(:,1),data(:,2),data(:,1),data(:,3),data(:,1),data(:,4),'-');

% labels
xlabel('Time / [dim. unit]');
ylabel('Displacement / [dim. unit]');

% legend
legend('Atom 1','Atom 2','Atom 3');

% axis limits
xlim([0,50]);
