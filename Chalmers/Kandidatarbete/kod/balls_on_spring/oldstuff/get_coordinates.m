function [ x, y, z, vx, vy, vz ] = get_coordinates( Q, dim )
%Gets the position and velocity vectors from a general coordinate Q.
%
%The format Q must have is:
%Q=[x1, x1', (y1, y1', z1, z1'), ... xN, xN', (yN, yN', zN, zN')]
%where ' denotes time derivative.


if dim<1 || dim>3
    disp('Error: dimention must be 1, 2 or 3. ')
    return
end


step=2*dim; %step to get to the next of the desiered coordinate

x =[];y =[];z =[];
vx=[];vy=[];vz=[];
if dim>=1;
    x =Q(:, 1:step:end);
    vx=Q(:, 2:step:end);
end

if dim>=2;
    y =Q(:, 3:step:end);
    vy=Q(:, 4:step:end);
end

if dim==3
    z =Q(:, 5:step:end);
    vz=Q(:, 6:step:end);
end

end

