function [ R_cm ] = get_CM_v2( Q, m, dim )
%Locates center of mass.
%result is a Lx1xdim 3D array.


R=get_coordinates_v3(Q,dim);
L=size(Q,1);

m_matrix=repmat(m, L,1,dim);
M=sum(m);%total mass

R_cm= sum(R.*m_matrix, 2)/M;



end

