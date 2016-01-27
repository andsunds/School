function [ Dq ] = dq_SDE(q, m, k, dampening, R0, N, dim )

q=reshape(q, 2,[]);%makes a matrix from the vector

r=reshape(q(1,:), dim,[])';


f=springforce_SDE( r,  k, R0 );


a=eye(N,N-1)-[zeros(1,N-1); eye(N-1)];%perhaps make this a spares matrix

F=a*f;

A=F./repmat(m, dim, 1).';

dampening=repmat(dampening, 1, dim);

Dq=[q(2,:);
     reshape(A.', 1, [])-dampening.*q(2,:)];

Dq=reshape(Dq, [], 1);

