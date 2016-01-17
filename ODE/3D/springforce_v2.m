function [ Fab ] = springforce_v2( r, k, R0 )
%The spring force (vector) acting on particle a.
% r = matrix of pos. vector of particle (row vector)
% k  = spring constant
% R0 = equilibrium distance
%This script can handle higer dimentions.

%dr=rb-ra;%vector

dr=diff(r);

Fab=dr.*k.*(1-R0./sqrt(sum(dr.^2, 2)));

end

