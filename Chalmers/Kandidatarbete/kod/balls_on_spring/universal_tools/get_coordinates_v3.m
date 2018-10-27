function [ R, V ] = get_coordinates_v3( Q, dim )
%Gets the position and velocity vectors from a general coordinate Q.
%Results are given as 3D arrays where all the particles x-components are in
%positions (:,:,1), y-components in (:,:,2), and z-components in (:,:,3).
%       time^                  particle^          dimention (in space)^
%
%The format Q must have is:
%Q=[x1, x1', (y1, y1', z1, z1'), xN, xN', (yN, yN', zN, zN'); @t=0
%   ........................................................; 
%   x1, x1', (y1, y1', z1, z1'), xN, xN', (yN, yN', zN, zN')] @t="end"
%where ' denotes time derivative.


if dim<1 || dim>3
    fprintf('Error: dimention must be 1, 2 or 3.\nGot dim = %d.\n\n',dim);
    R=[];V=[];
    return
end

step=2*dim; %step to get to the next of the desiered coordinate

%initilizing
R= zeros([size(Q(:, 1:step:end)), dim]);
V= zeros([size(Q(:, 2:step:end)), dim]);


for i=1:dim
    R(:,:,i)= Q(:, (2*i-1):step:end);
    V(:,:,i)= Q(:,   (2*i):step:end);
end


end

