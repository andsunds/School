function [ Dq ] = dq_two_balls(q, m1, m2, k, n, R0, Br )


q=reshape(q, 2,[]);%makes a matrix from the vector


r1=q(1,1:2);
r2=q(1,3:4);
%r1=[q(1), q(7)];
%r2=[q(4), q(10)];

F12=springforce_v1( r1, r2,  k, R0 );

A=[F12/m1, -F12/m2];

%Dq=[q(2), q(3), A(1), q(5),q(6), A(2), q(8),q(9), A(3), q(11),q(12), A(4)].';

 Dq=[q(2,:)+Br*randn(1,4);
     A-n*q(2,:)];

Dq=reshape(Dq, [], 1);

