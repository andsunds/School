function [  ] = play_movie_v3(T,Q,m,dim,framesize, playbackspeed, followCM)
%Plays a movie sequence from a simulated event
% T:   vector with all the time instances
% Q:   position and velocity vector
% m:   vector containing all the particles masses
% dim: dimentions, must be 2 or 3
% followCM: 1 or 0
%
%The format Q must have is:
%Q=[x1, x1', y1, y1', (z1, z1'), ... xN, xN', yN, yN', (zN, zN')]
%where ' denotes time derivative.

if followCM~=0
    followCM=1;
end


%    v-start with a pause
dt=[.2; diff(T)./playbackspeed];%determines pause time between frames

%different colors for each particle
Colors=hsv(length(m));


R   =get_coordinates_v3(Q,dim);
R_cm=get_CM_v2(Q,m,dim);

switch dim
    case 1 %improve for 1D
        t=repmat(T,1,length(m));
        p=   plot(t(1,:), R(1,:,1),'-'); hold on
        s=scatter(t(1,:), R(1,:,1), 81, Colors); 
        h=   plot(t(1,:),R_cm(1,:,1));
        axis equal; grid on
        xlabel('$t$', 'Interpreter', 'Latex', 'FontSize', 20, 'Color', 'k');
        ylabel('$x$', 'Interpreter', 'Latex', 'FontSize', 20, 'Color', 'k');
        set(gca,'FontSize',15);
        for i=1:(length(T)-1)
            pause(dt(i))  
            set(p, 'XData',t(i,:),'YData',R(i,:,1));
            %plot(t(1:i,:), R(1:i,:,1),'-k');
            set(s, 'XData',t(i,:)  ,'YData',R(i,:,1));
            set(h, 'XData',T(i,:),'YData',R_cm(i,1,1));
            axis([T(i)*followCM-framesize, T(i)+framesize,...
                R_cm(i,1,1)*followCM-framesize, R_cm(i,1,1)*followCM+framesize]);
        end%end for
        
    case 2
        p=   plot(R(1,:,1),R(1,:,2),'-k'); hold on
        s=scatter(R(1,:,1),R(1,:,2), 81, Colors); 
        h=plot(R_cm(1,1,1),R_cm(1,1,2), 'xr');
        axis equal; grid on
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 20, 'Color', 'k');
        ylabel('$y$', 'Interpreter', 'Latex', 'FontSize', 20, 'Color', 'k');
        set(gca,'FontSize',15);
        for i=1:(length(T)-1)
            pause(dt(i))  
            set(p, 'XData',R(i,:,1),'YData',R(i,:,2));
            set(s, 'XData',R(i,:,1),'YData',R(i,:,2));
            set(h, 'XData',R_cm(i,1,1),'YData',R_cm(i,1,2));
            axis([R_cm(i,1,1)*followCM-framesize, R_cm(i,1,1)*followCM+framesize,...
                  R_cm(i,1,2)*followCM-framesize, R_cm(i,1,2)*followCM+framesize]);
        end%end for
        
    case 3
        p=   plot3(R(1,:,1),R(1,:,2),R(1,:,3),'-k'); hold on
        s=scatter3(R(1,:,1),R(1,:,2),R(1,:,3), 81, Colors); 
        h=plot3(R_cm(1,1,1),R_cm(1,1,2),R_cm(1,1,3), 'xr');
        axis equal; grid on
        xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 20, 'Color', 'k');
        ylabel('$y$', 'Interpreter', 'Latex', 'FontSize', 20, 'Color', 'k');
        zlabel('$z$', 'Interpreter', 'Latex', 'FontSize', 20, 'Color', 'k');
        set(gca,'FontSize',15);
        for i=1:(length(T)-1)
            pause(dt(i))  
            set(p, 'XData',R(i,:,1),'YData',R(i,:,2),'ZData',R(i,:,3));
            set(s, 'XData',R(i,:,1),'YData',R(i,:,2),'ZData',R(i,:,3));
            set(h, 'XData',R_cm(i,1,1),'YData',R_cm(i,1,2),'ZData',R_cm(i,1,3));
            axis([R_cm(i,1,1)*followCM-framesize, R_cm(i,1,1)*followCM+framesize,...
                  R_cm(i,1,2)*followCM-framesize, R_cm(i,1,2)*followCM+framesize,...
                  R_cm(i,1,3)*followCM-framesize, R_cm(i,1,3)*followCM+framesize]);
        end%end for
end%end switch

hold off

end


