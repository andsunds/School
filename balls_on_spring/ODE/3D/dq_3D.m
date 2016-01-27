function [ Dq ] = dq_3D(q, m, k, n, R0, N, dim )


q=reshape(q, 2,[]);%makes a matrix from the vector

r=reshape(q(1,:), dim,[])';


f=springforce_v2( r,  k, R0 );


a=eye(N,N-1)-[zeros(1,N-1); eye(N-1)];

F=a*f;

A=F./repmat(m, dim, 1).';

Dq=[q(2,:);
     reshape(A.', 1, [])-n*q(2,:)];

Dq=reshape(Dq, [], 1);

