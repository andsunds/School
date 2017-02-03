function [ T ] = spTranferMatrix(J, alpha, beta)
%Creates a sparse matrix that corresponds to one time step.
%   The matrix is tri-diagonal and has non-zero elements in the
%   anti-diagonal corners. 

D0=1-beta;
D1=(beta-alpha)/2;

%initializing the sparese transfermatrix
T=spalloc(J,J,3*J);
B=repmat([-D1, D0, +D1],J,1);
T=spdiags(B,-1:1,T);
T(1,end)=-D1; %Periodic BC
T(end,1)=+D1; %Periodic BC
%spy(T) %debug

end

