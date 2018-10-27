function [ celsius ] = V2C( volt, order )


files=dir('Data/kalibrering_v2_*C.txt');
l=length(files);

T=zeros(l,1);
V=zeros(l,5);

for i=1:l
    data=load(files(i).name);
    T(i)=str2double(regexp(files(i).name(15:end), '\d+', 'match'));
    V(i,:)=mean(reshape(data,5,[]),2).';
end

V=V(:,order);

k=zeros(2,5);
for i=1:size(k,2)
    k(:,i)=[V(:,i),ones(3,1)]\T;
end


celsius=zeros(size(volt));
for i=1:size(volt,2)
    celsius(:,i)=k(1,i)*volt(:,i)+k(2,i);
end






end

