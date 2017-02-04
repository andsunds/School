function [ T ] = spTranferMatrix(J, alpha, beta)
%Creates a sparse matrix that corresponds to one time step.
%   The matrix is tri-diagonal and has non-zero elements in the
%   anti-diagonal corners. 

%The values of the diagonal (+-1) elements
D0=1-beta;
Dp1=(beta-alpha)/2; 
Dm1=(beta+alpha)/2;

%initializing the sparese transfermatrix
T=spalloc(J,J,3*J);
B=repmat([Dm1, D0, Dp1],J,1);
T=spdiags(B,-1:1,T);
T(1,end)=Dm1; %Periodic BC
T(end,1)=Dp1; %Periodic BC

end