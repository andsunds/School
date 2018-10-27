function [ I ] = create_indecis( N )
%Tar fram index sorterade så att de övre diagonalerna hamnar på respektive
%rad. 
%OBS: bör bara användas för en uppåt triangulär matris.
%
%Exempel:
%[ * ~ + x » ]          [ * * * * * ]
%| 0 * ~ + x |          | ~ ~ ~ ~ 0 |
%| 0 0 * ~ + |   --->   | + + + 0 0 |
%[ 0 0 0 * ~ |          | x x 0 0 0 |
%[ 0 0 0 0 * ]          [ » 0 0 0 0 ]
%

%Skapar en matris där varje rad inehåller diagonalens index.
I=repmat(1:(N+1):N^2,N,1);
%Här brys upprepningarna längs raderna, nu får varje rad (nästan) indexena 
%för varje respektive avdiagonal
I=I+N*repmat((0:(N-1)).', 1,N);
%Men det blir lite för stora värden i ändarna på raderna. 
%Det åtgärdas här:
I(I>N^2)=I(I>N^2) - N^2;

end
