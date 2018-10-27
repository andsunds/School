function [ TN ] = koordinatbyte( XY )


n=minsta_kvadrat( XY(:,1), XY(:,2) ); %normalen

n=n/norm(n);%normera normalen
t=[n(2); -n(1)];%tangentvektor

TN=XY*[t, n];%koordinatbyte



end

