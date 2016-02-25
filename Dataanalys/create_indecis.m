function [ I ] = create_indecis( N )
%Creates a matrix of indecis for the upper diagonals.

I=repmat(1:(N+1):N^2,N,1);

I=I+N*repmat((0:(N-1)).', 1,N);

I(I>N^2)=I(I>N^2) - N^2;

end
