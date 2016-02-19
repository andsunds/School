clc, clear all
A=importdata('confined_270304-6-28-min.mat');

%%
C=cell(size(A,3),1);
XY=cell(size(A,3),1);
for i=1:size(A,3)
    clearvars rad kol gradX gradY
    [rad,kol] = find(255==A(:,:,i));
    XY{i}=[rad,kol];
    grad= diff(XY{i});
    C{i}=grad;
    %= abs(diff(i)-diff(i-1))
end

%%
clf
res=5
[rad,kol] = find(255==A(:,:,1));
X=length(kol);
grid=floor(X/res);
for i=1:res
    x(i)=grid*i;
    y(i)=rad(x(i))
end
p = polyfit(x,y,res-1);

x1 = linspace(min(kol),x(res));
y1 = polyval(p,x1);
plot(x1,y1)


