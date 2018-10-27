function [  ] = play_movie(T,Q,m,dim,framesize, playbackspeed )
% T:   vector with all the time instances
% Q:   position and velocity vector
% m:   vector containing all the particles masses
% dim: dimentions, must be 2 or 3
%
%The format Q must have is:
%Q=[x1, x1', y1, y1', (z1, z1'), ... xN, xN', yN, yN', (zN, zN')]
%where ' denotes time derivative.



dt=diff(T)./playbackspeed;
L=length(T);
N=length(m);
M=sum(m);

%different colors for each particle
%c=logspace(-2, 0, N).';
%C=[c, circshift(c, [1 0]), circshift(c, [2 0])];
Colors=hsv(N);


%position of particles
step=2*dim;
x=Q(:, 1:step:end);
y=Q(:, 3:step:end);


%center of mass pos.
m_matrix=repmat(m, L,1);
x_cm=sum(m_matrix.*x, 2)/M;
y_cm=sum(m_matrix.*y, 2)/M;


if dim==3
    z=Q(:, 5:step:end);
    z_cm=sum(m_matrix.*z, 2)/M;
    %'bonds'
    p=plot3(x(1,:),y(1,:),z(1,:),'-k'); hold on
    %individual particles
    s=scatter3(x(1,:),y(1,:),z(1,:), 81, Colors); hold on
    %center of mass
    h=plot3(x_cm(1),y_cm(1), z_cm(1),'xr');
else
    p=plot(x(1,:),y(1,:),'-k'); hold on
    s=scatter(x(1,:),y(1,:), 81, Colors); hold on
    h=plot(x_cm(1),y_cm(1), 'xr');
end%end if


axis equal
grid on


if dim == 3
    for i=1:L-1
        set(p, 'XData',x(i,:),'YData',y(i,:),'ZData',z(i,:));
        set(s, 'XData',x(i,:),'YData',y(i,:),'ZData',z(i,:));
        set(h, 'XData',x_cm(i),'YData',y_cm(i),'ZData',z_cm(i));
        axis([x_cm(i)-framesize, x_cm(i)+framesize,...
              y_cm(i)-framesize, y_cm(i)+framesize, ...
              z_cm(i)-framesize, z_cm(i)+framesize]);
        pause(dt(i)) %this only has length L-1
    end%end for
    %need special case for the last frame
    i=i+1;
    set(p, 'XData',x(i,:),'YData',y(i,:),'ZData',z(i,:));
    set(s, 'XData',x(i,:),'YData',y(i,:),'ZData',z(i,:));
    set(h, 'XData',x_cm(i),'YData',y_cm(i),'ZData',z_cm(i));
    axis([x_cm(i)-framesize, x_cm(i)+framesize,...
          y_cm(i)-framesize, y_cm(i)+framesize, ...
          z_cm(i)-framesize, z_cm(i)+framesize]);
    
else
    for i=1:L-1
        set(p, 'XData',x(i,:),'YData',y(i,:));
        set(s, 'XData',x(i,:),'YData',y(i,:));
        set(h, 'XData',x_cm(i),'YData',y_cm(i));
        axis([x_cm(i)-framesize, x_cm(i)+framesize,...
              y_cm(i)-framesize, y_cm(i)+framesize]);
        pause(dt(i))
    end%end for
    %need special case for the last frame
    i=i+1;
    set(p, 'XData',x(i,:),'YData',y(i,:));
    set(s, 'XData',x(i,:),'YData',y(i,:));
    set(h, 'XData',x_cm(i),'YData',y_cm(i));
    axis([x_cm(i)-framesize, x_cm(i)+framesize,...
          y_cm(i)-framesize, y_cm(i)+framesize]);
end%end if

end%end function

