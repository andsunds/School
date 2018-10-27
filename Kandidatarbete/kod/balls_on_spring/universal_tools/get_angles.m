function [angles] = get_angles(Q,dim)


R=get_coordinates_v3(Q,dim);


S=size(R);
R=cat(3,R,zeros([S(1:2),3-dim]));

dR=diff(R,1,2);

%product of norms of two consecutive pos. vectors
dR_norm_prod=(sqrt(sum(dR(:,1:end-1,:).^2, 3)).*sqrt(sum(dR(:,2:end,:).^2, 3)));


xprod=cross(dR(:,1:end-1,:),dR(:,2:end,:),3);
sin_angle=sqrt(sum(xprod.^2,3))./dR_norm_prod;

dotprod=dot(dR(:,1:end-1,:),dR(:,2:end,:),3);
cos_angle=dotprod./dR_norm_prod;

angles=atan(sin_angle./cos_angle);
angles(angles<0)=angles(angles<0)+pi;

grid on


