function [ I ] = create_indecis( N )

I=repmat(1:(N+1):N^2,N,1);

I=I+repmat((0:(N-1)).'*N, 1,N);

I(I>N^2)=I(I>N^2) - N^2;

end
