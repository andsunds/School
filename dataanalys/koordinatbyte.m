function [ TN ] = koordinatbyte( XY )


n=minsta_kvadrat( XY(:,1), XY(:,2) ); %normalen

n=n/norm(n);
t=[n(2); -n(1)];

TN=XY*[t, n];



end

