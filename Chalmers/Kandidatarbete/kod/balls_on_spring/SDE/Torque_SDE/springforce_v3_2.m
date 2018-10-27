function [ F ] = springforce_v3_2( r,  k, R0, s, theta0, dim )
%The total force (vector) acting on every particle. 
% r      = matrix of pos. vector of particle (row vector)
% k      = radial spring constants, (N-1)x1 vector
% R0     = equilibrium distance, (N-1)x1 vector 
% s      = axial spring constants, (N-2)x1 vector
% theta0 = equilibrium angle, (N-2)x1 vector 
%This script can handle higer dimentions.
%_NOT_ backwards compatible!!!

N=size(r,1);

r=[r,zeros(N,3-dim)];%pading to ensure 3D


dr=diff(r, 1, 1);


%radial force
%f_radial=dr.*k.*( 1-R0./sqrt(sum(dr.^2, 2)) );
f_radial=dr.*k.*( 1-R0./row_norm(dr, R0) - .5*R0.^3./row_norm(dr, R0).^3 - row_norm(dr, R0)./R0 );

%this matrix distributes the tensions i all springs as forces acting on
%each particle
a=eye(N,N-1)-[zeros(1,N-1); eye(N-1)];
F_radial=a*f_radial;



%axial force

%product of norms of two consecutive pos. vectors
dr_norm_prod= row_norm(dr(1:end-1,:), R0) .* row_norm(dr(1:end-1,:), R0);

%sine
xprod=cross(dr(1:end-1,:),dr(2:end,:),2);
sin_angle=row_norm(xprod, R0.^2)./dr_norm_prod;
%cosine
dprod=dot(dr(1:end-1,:),dr(2:end,:),2);
cos_angle=dprod./dr_norm_prod;
%tan
angle=atan(sin_angle./cos_angle);%no risk of invalid angels due to |sin_angle|>1
angle(angle<0)=angle(angle<0)+pi;


%         v-normalize angle has been calculated
Mo=(xprod./row_norm(xprod, R0.^2)).*tan((angle-theta0)/1).*s;
%                                         v--concatenate with same force on
%                                         v  the next particle
f_axial=reshape([cross(dr(1:end-1,:),Mo,2), cross(dr(2:end,:),Mo,2)].', 3,[]).';
f_axial=reshape(f_axial.', 3,[]).';%this makes f_axial alternate each row as above


%creates maxtrix to distribute torque forces cf. 'a'
b=[repmat([1 -1 -1 1], N,1), zeros(N, size(f_axial,1)-1)];
for i=0:N-1
    b(1+i,:)=circshift(b(i+1,:), i*2,2);
end
b=b(:,4:end);

F_axial=b*f_axial;


F=F_axial(:,1:dim)+F_radial(:,1:dim);%removes padding (if needed)

end

function [Norm] = row_norm(A, characteristic_lenght)
Norm=sqrt(sum(A.^2, 2));
%If we get 0, 
%then set it to something very small compared to the characteristic length
Norm(Norm==0)=mean(characteristic_lenght)*1e-9;
end