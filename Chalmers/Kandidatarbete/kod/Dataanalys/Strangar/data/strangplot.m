%% strängplot
%kod för att generera plottdata till gnuplot
clc;clf;clearvars


filnamn1=cell(1,4);
filnamn1{1}='confined_28min.mat'; 
filnamn1{2}='confined_32min.mat';
filnamn1{3}='nonconfined_5min.mat';
filnamn1{4}='nonconfined_167min.mat';

filnamn2=cell(1,4);
filnamn2{1}='confined_28min_polynom.mat'; 
filnamn2{2}='confined_32min_polynom.mat';
filnamn2{3}='nonconfined_5min_polynom.mat';
filnamn2{4}='nonconfined_167min_polynom.mat';

fil=4;
load([filnamn1{fil}]);
load([filnamn2{fil}], 'px', 'py');

i=50;
%polynom
I=1000;
s=linspace(0,1,I);
x=polyval(px(i,:),s);
y=polyval(py(i,:),s);
%spline
t=1:N_points(i);
xy = [coordinates(i,t,1); coordinates(i,t,2)];
sp=csape(t,xy,'spline');
XY=ppval(linspace(0,365,I), sp);

%plottning
plot(coordinates(i,:,1),coordinates(i,:,2), '-')
hold on
plot(x,y)
plot(XY(1,:), XY(2,:))

axis equal
axis([16,17,9.2,9.7]*1e-6)

%%
clc
plotdata=NaN(1000, 6);

plotdata(1:length(coordinates(i,:,1)), 1:2) = [coordinates(i,:,1); coordinates(i,:,2)].';
plotdata(:,3:6)=[x;y;XY].';

save('strang_anpassning.tsv', 'plotdata', '-ascii')
