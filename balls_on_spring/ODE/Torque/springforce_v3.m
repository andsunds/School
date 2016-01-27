function [ F ] = springforce_v3( r,  k, R0, s, theta0, dim )
%The total force (vector) acting on every particle. 
% r      = matrix of pos. vector of particle (row vector)
% k      = radial spring constants, (N-1)x1 vector
% R0     = equilibrium distance, (N-1)x1 vector 
% s      = axial spring constants, (N-2)x1 vector
% theta0 = equilibrium angle, (N-2)x1 vector 
%This script can handle higer dimentions.
%_NOT_ backwards compatible!!!

N=size(r,1);

r=[r,zeros(N,3-dim)];%pading in case of lowe dimention

dr=diff(r, 1, 1);


%radial force
%k=reshape(k, [],1);
%R0=reshape(R0, [],1);
f_radial=dr.*k.*( 1-R0./sqrt(sum(dr.^2, 2)) );
%f_radial=dr.*k.*( 1-R0./row_norm(dr, R0) );

%this matrix distributes the tensions i all springs as forces acting on
%each particle
a=eye(N,N-1)-[zeros(1,N-1); eye(N-1)];
F_radial=a*f_radial;


%axial force
%s=reshape(s, [],1);
%theta0=reshape(theta0, [],1);
xprod=cross(dr(1:end-1,:),dr(2:end,:),2);

%sin_angle=sqrt(sum(xprod.^2,2))./(sqrt(sum(dr(1:end-1,:).^2, 2)).*sqrt(sum(dr(2:end,:).^2, 2)));
sin_angle=row_norm(xprod, R0.^2)./( row_norm(dr(1:end-1,:), R0).*row_norm(dr(1:end-1,:), R0) );
%angle=asin(sin_angle);

%cos_angle=dot(dr(1:end-1,:),dr(2:end,:),2)./(sqrt(sum(dr(1:end-1,:).^2, 2)).*sqrt(sum(dr(2:end,:).^2, 2)));
cos_angle=dot(dr(1:end-1,:),dr(2:end,:),2)./( row_norm(dr(1:end-1,:), R0).*row_norm(dr(1:end-1,:), R0) );


%angle=asin(sin_angle);
%angle(sin_angle>=1)=acos(cos_angle);
%angle(sin_angle<=1)=acos(cos_angle);
%angle(cos_angle<0)=pi-angle(cos_angle<0);
angle=atan(sin_angle./cos_angle);
angle(angle<0)=angle(angle<0)+pi;




xprod=xprod./row_norm(xprod, R0.^2);%normalize AFTER angle has been calculated

Mo=xprod.*(angle-theta0).*s;
f_axial=reshape([cross(dr(1:end-1,:),Mo,2),cross(dr(2:end,:),Mo,2)].', dim,[]).';
f_axial=reshape(f_axial.', dim,[]).';


%creates maxtrix to distribute torque forces cf. 'a'
b=[repmat([1 -1 -1 1], N,1), zeros(N, size(f_axial,1)-1)];
for i=0:N-1
    b(1+i,:)=circshift(b(i+1,:), i*2,2);
end
b=b(:,4:end);

F_axial=b*f_axial;

%pause
F=F_axial(:,1:dim)+F_radial(:,1:dim);

%F=a*f_radial + b*f_axial;



end


function [Norm] = row_norm(A, car_len)
%gets the norm along the rows
Norm=sqrt(sum(A.^2, 2));
%If we get 0, 
%then set it to something very small
%compared to the characteristic length R0:
Norm(Norm==0)=mean(car_len)*1e-9;
end

% function [normalized] = normalize(A, skip_0)
% %only for use in this script
% % A is a matrix containing row vectors
% nor=row_norm(A);
% 
% 
% 
% 
% normalized=A./nor;
% end










