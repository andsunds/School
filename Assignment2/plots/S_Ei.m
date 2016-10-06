function [ Sn ] = S_Ei( n, x, a )


l=length(x);

A=repmat([1,cumprod(a-(0:(n-1)))].', 1,l);

SX=repmat(x,n+1,1).^repmat((0:n).' ,1 , l);

%size(SX)
%size(A)

Y=A./SX;
%size(Y)


Sn= (exp(-x).*x.^a).*sum(Y,1);


end


