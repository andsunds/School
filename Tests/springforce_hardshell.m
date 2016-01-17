function [ Fab ] = springforce_v4( r, k, R0 )
%The spring force (vector) acting on particle a.
% r = matrix of pos. vector of particle (row vector)
% k  = spring constant
% R0 = equilibrium distance
%This script can handle higer dimentions.

%dr=rb-ra;%vector

dr=diff(r);

spring=dr.*k.*(1-R0./sqrt(sum(dr.^2, 2)));

hard_shell=-dr.*(R0/5-sqrt(sum(dr.^2, 2))).^(-4);

Fab=spring+hard_shell;

end

