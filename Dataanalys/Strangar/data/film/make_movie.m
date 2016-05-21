%%
clc;
clear all;
clf

filnamn=cell(1,4);
filnamn{1}='confined_28min_polynom.mat'; 
filnamn{2}='confined_32min_polynom.mat';
filnamn{3}='nonconfined_5min_polynom.mat';
filnamn{4}='nonconfined_167min_polynom.mat';

fil=4;

load(['../', filnamn{fil}]);

i=1;

S=linspace(0,1,1000);

x=polyval(px(i,:), S);
y=polyval(py(i,:), S);

plot(polyval(PX_mean, S),polyval(PY_mean, S), 'k', 'linewidth', 4)
hold on


scalebar_x=[1.45,1.45,1.45, 2.45,2.45,2.45]*1e-5;
scalebar_y=[.45,.55,.5, .5,.45,.55]*1e-5;
plot(scalebar_x,scalebar_y,'k')
text(1.9e-5,.57e-5,'10 µm', 'fontsize',15)
axis equal
axis([0,2.5, .4,2.2]*1e-5)

set(gca,'XTick',[],'yTick',[]);
%set(gca,'yTick',[]);


h=plot(x,y, 'b', 'linewidth', 2);
axis([0,2.5, .4,2.2]*1e-5)
axis equal

v=VideoWriter('filamen4.avi');
v.FrameRate = 10;
open(v)

for i=1:100
    set(h, 'XData', polyval(px(i,:), S),'yData', polyval(py(i,:), S))
    axis([0,2.5, .4,2.2]*1e-5)
    pause(.1)
    F = getframe;
    writeVideo(v,F)
end
close(v)