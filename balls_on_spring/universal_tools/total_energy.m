function [ E ] = total_energy( Q, m, k, dim )
%Calculates total energy in a system of particles with masses given by the
%vector 'm' with springs in between with sprin constants given by the
%vector 
%
%The format Q must have is:
%Q=[x1, x1', (y1, y1', z1, z1'), xN, xN', (yN, yN', zN, zN'); @t=0
%   ........................................................; 
%   x1, x1', (y1, y1', z1, z1'), xN, xN', (yN, yN', zN, zN')] @t="end"
%where ' denotes time derivative.


[R, V] = get_coordinates_v3( Q, dim );
L=size(Q,1);

%                  v-make sure m is a row vector
m_matrix=repmat(reshape(m,1,[]), L, 1);

E_kinetic=0.5*sum(sum(V.^2, 3).*m_matrix, 2);
%sum over:   all dimentions-^             ^-all particles


%                  v-make sure k is a row vector
k_matrix=repmat(reshape(k,1,[]), L, 1);

%         v-first order diff
dR=diff(R,1,2);
%           ^-over the individual particles

E_spring=0.5*sum(sum(dR.^2, 3).*k_matrix, 2);
%sum over:   all dimentions-^             ^-all particles

E=E_spring+E_kinetic;

end

