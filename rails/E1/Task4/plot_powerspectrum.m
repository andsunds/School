% plot the powerspectrum
% Created by Martin Gren 2014-10-25.

% load the data file
data = importdata('powerspectrum.dat');

%plot
figure;
plot(data(:,1),data(:,2),'-');

% labels
xlabel('Frequency','Interpreter','LaTeX');
ylabel('Powers pectrum','Interpreter','LaTeX');



