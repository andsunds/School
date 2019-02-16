% plot the powerspectrum
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('function.dat');

%plot
figure;
plot(data(:,1),data(:,2),'-');

% labels
xlabel('$t$','Interpreter','LaTeX');
ylabel('$h(t)$','Interpreter','LaTeX');



