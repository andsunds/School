clc, clear all
A=importdata('confined_270304-6-28-min.mat');

%% Gradient 
C=cell(size(A,3),1);
XY=cell(size(A,3),1);
for i=1:1
    clearvars rad kol gradX gradY
    [rad,kol] = find(255==A(:,:,i));
    XY{i}=[rad,kol];
    grad= diff(XY{i});
    C{i}=grad;
%     X=abs(grad(:,1))>0;
%     lut=find(1==X);
%     for j=2:length(lut)
%         dx(j)=grad(lut(j),1);%/sum(grad(lut(j-1):1:lut(j),2));
%         dy(j)=sum(grad(lut(j-1):1:lut(j),2));
%     end
%     
end

%% polynomanpassning
clf
for j=1:size(A,3)
[rad,kol] = find(255==A(:,:,j));    %hämtar filament nummer j
res=length(kol);                    %upplösning i 'x-led'
clearvars y1
for i=1:100
grad=i;                             %polynomanpassning av grad i
p = polyfit(kol,rad,grad);          
x1 = linspace(min(kol),max(kol),res);
y1(i,:) = polyval(p,x1);
fel(i)=sum(abs(y1(i,:)'-rad));      %bästa anpassningen = minsta felet
end
[a,bragrad]=find(min(fel)==fel);    
plot(x1,y1(bragrad,:),kol,rad)    
legend(sprintf('gradP=%0.2f',bragrad))
pause(.2)                           
end

%% 






