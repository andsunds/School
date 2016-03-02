function [ tangent, normal ] = tangent_normal( polynom_x , polynom_y, punkter )
%Tar polynomen för x och y samt i vilka punkter som tangent och normal ska
%undersökas. 
%
%Ger ut normerad tangent- och normalvektor i dessa punkter. 
%Formatet på utskrifterna är en 3D-matris med varje tidssteg i den tredje
%dimensionen. 
%
%x-komp finns i 1: raden i 1:a dimensionen, y i 2:a raden i
%1:a dimensionen. 
%
%Varje punkts resp. värden finns längs med 2:a dimensinen.


N=size(polynom_x, 1);

tangent=zeros(2, length(punkter) , N);%init.

for i=1:N
    %Tangenten ges av derivatan:
    dx=polyder(polynom_x(i,:));%derivera x
    dy=polyder(polynom_y(i,:));%derivera y
    
    %Tangentvektor i de specifika pkt.
    tmp=[ polyval(dx, punkter); polyval(dy, punkter) ];
    %Varje tidsögobicks tangenter.
    tangent(:,:,i)=tmp./repmat(sqrt(sum(tmp.^2,1)),2,1);%normering
end
%Normalen fås enkelt från tangenten.
normal=[tangent(2,:,:); -tangent(1,:,:)];


end

