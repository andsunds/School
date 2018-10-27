function [ Dq ] = dq_torque(q, m, k, R0, s, theta0, damp, dim )
%Calculates the differential of a generalized position vector q.
% m   = vector containing all the particles masses
% k   = vector with all radial spring constants
% s   = vector with all axial spring constants
% damp= dampening koefficient



q=reshape(q, 2,[]);%makes a matrix from the vector

r=reshape(q(1,:), dim,[]).';


F=springforce_v3_2( r,  k, R0, s, theta0, dim);



A=F./repmat(m, dim, 1).';

Dq=[q(2,:);
     reshape(A.', 1, [])-damp*q(2,:)];

Dq=reshape(Dq, [], 1);

