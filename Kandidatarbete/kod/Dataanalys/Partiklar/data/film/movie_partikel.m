%%
clc;
clear all;
clf

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv'; 
filnamn{2}='logphasecells.csv';


fil=1;

data=load(['../', filnamn{fil}]);
addpath('../../');

C=separera(data);

part1=75;
part2=167;
part3=205;

XY1=C{part1}(:,2:3);
I1=mean(C{part1}(:,4));
XY2=C{part2}(:,2:3);
I2=mean(C{part2}(:,4));
XY3=C{part3}(:,2:3);
I3=mean(C{part3}(:,4));

%max(abs(XY1))

i=1;
hold on
k=50;
h1=plot(XY1(i,2),XY1(i,1), 'b.', 'markersize',I1^(1/3)*k);
h2=plot(XY2(i,2),XY2(i,1), 'r.', 'markersize',I2^(1/3)*k);
h3=plot(XY3(i,2),XY3(i,1), 'g.', 'markersize',I3^(1/3)*k);
axis([-1.8,1.8, -1,1]*1e-7)
axis equal
axis([-1.8,1.8, -1,1]*1e-7)





scalebar_x=[.7,.7,.7, 1.7,1.7,1.7]*1e-7;
scalebar_y=-[.85,.95,.9, .9,.85,.95]*1e-7;
plot(scalebar_x,scalebar_y,'k')
text(1.05e-7,-0.83e-7,'100 nm', 'fontsize',15)
set(gca,'XTick',[],'yTick',[]);

%axis equal
%axis([0,2.5, .4,2.2]*1e-5)

%% Skapar en film

v=VideoWriter('partiklar_ed.avi');
v.FrameRate = 7;
open(v)

for i=1:70
    set(h1, 'XData', XY1(i,2),'yData',XY1(i,1) )
    set(h2, 'XData', XY2(i,2),'yData',XY2(i,1) )
    set(h3, 'XData', XY3(i,2),'yData',XY3(i,1) )
    axis([-1.8,1.8, -1,1]*1e-7)
    pause(1/v.FrameRate)
    F = getframe;
    writeVideo(v,F)
end
close(v)

%% För att ta reda på vilka index man ska välja
clc;clf;clearvars
load('../kompleterande_data.mat')

%plot(intensitet{1}, '.')
plot(sqrt(std_n{1}.^2+std_t{1}.^2), '.')

