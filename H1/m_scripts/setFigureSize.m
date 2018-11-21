function [ fig ] = setFigureSize( fig )
%figureSizePaper1 
fig.Units = 'points';
W = 600; 
H = 300;

fig.WindowStyle = 'normal'; % undock 
fig.Position(3:4) = [W H];

end
