function [ ret ] = RHS_ignition( TV, nuT, nun )
% integrating from s=0 to s=1, s=1-r^2/a^2
I=integral(@(s) sigma_v(TV, nuT, s).*s.^(2*nun), 0,1,'RelTol',1e-16);
%Calculated K_\alpha * P * \tau_E / K_L
ret=(1+2*nun)*(1+nuT)^2/((1+nun+nuT)*(1+nun)^2).*TV.^2./I;
end


function [ret] = sigma_v(TV, nuT, s)

%coefficients
alpha=-0.2935; 
am1=-21.38; a0=-25.2; a1=-7.101e-2;
a2=1.938e-4; a3=4.925e-6;a4=-3.984e-8;

%T(s)
Ts=(1+nuT)*TV*(s.^nuT);

% <\sigma v>(T)
ret=1e16*exp(am1*Ts.^alpha + a0 + a1*Ts + a2*Ts.^2 + ...
             a3*Ts.^3 + a4*Ts.^4);
end
