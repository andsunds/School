function [ I ] = find_endpoint(picture, rad, kol )





L=length(rad);
C=zeros(L,1);
%B=zeros(5,5);

for i=1:length(rad)
    B=picture(rad(i)+(-2:2), kol(i)+(-2:2) );
    C(i)=length(find(B));
end

D=find(3==C);
reserv=find(4==C);

if length(D)==1 && kol(D)>kol(reserv(1))
    D=reserv;
    %disp('reserv, 1'), D(1)
elseif isempty(D)
    D=reserv;
    %disp('reserv, 0'), D(1)
end

%endpoint=reshape([kol(D(1)),rad(D(1))], 1,1,2);

I=D(1);

end

