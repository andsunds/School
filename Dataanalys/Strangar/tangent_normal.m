function [ tangent, normal ] = tangent_normal( polynom_x , polynom_y, punkter )
%Tar polynomen för x och y samt i vilka punkter som tangent och normal ska
%undersökas. 
%
%Ger ut normerad tangent- och normalvektor i dessa punkter. 


N=size(polynom_x, 1);

tangent=zeros(2, N);

for i=1:N
    %Tangenten ges av derivatan:
    dx=polyder(polynom_x(i,:));%derivera x
    dy=polyder(polynom_y(i,:));%derivera y
    
    %Tangentvektor i de specifika pkt.
    tmp=[polyval(dx, punkter);
         polyval(dy, punkter) ];
 
    tangent(i)=tmp./repmat(sqrt(sum(tmp.^2,1)),2,1);%normering
end

normal=[tangent(2,:); -tangent(1,:)];


end

