function [ r, v ] = get_coordinates_v2( Q, dim )
%Gets the position and velocity vectors from a general coordinate Q.
%Results are given as 3D arrays where all the particles x-components are in
%positions (:,:,1), y-components in (:,:,2), and z-components in (:,:,3).
%           ^all times                 ^all individual particles
%
%The format Q must have is:
%Q=[x1, x1', (y1, y1', z1, z1'), xN, xN', (yN, yN', zN, zN'); @t=0
%   ........................................................; 
%   x1, x1', (y1, y1', z1, z1'), xN, xN', (yN, yN', zN, zN')] @t="end"
%where ' denotes time derivative.


if dim<1 || dim>3
    disp('Error: dimention must be 1, 2 or 3. ')
    return
end


step=2*dim; %step to get to the next of the desiered coordinate

r= cat(3,Q(:, 1:step:end));
v= cat(3,Q(:, 2:step:end));


if dim>=2;
    r= cat(3,r,Q(:, 3:step:end));
    v= cat(3,v,Q(:, 4:step:end));
end

if dim==3
    r =cat(3,r,Q(:, 5:step:end));
    v =cat(3,v,Q(:, 6:step:end));
end


end

