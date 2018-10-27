function [ Fab ] = springforce_old( ra, rb, k, R0 )
%The spring force (vector) acting on particle a.
% ra = pos. vector of particle a (row vector)
% rb = pos. vector of particle b (row vector)
% k  = spring constant
% R0 = equilibrium distance
%This script can handle higer dimentions.

dr=rb-ra;%vector

Fab=dr.*k.*(1-R0./sqrt(sum(dr.^2, 2)));

end

