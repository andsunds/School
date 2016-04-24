%% Korrelation i r??relsen
clc;clf;clear all

load('filnamn.mat')

for fil=1:2 %to plot data for both cell types

data =load(filnamn{fil});
C = separera(data);
n=length(C);%antal partiklar

for i=1:n
    XY=C{i}(2:end-1,2:3);
    D1=diff(XY(1:end-1));
    D2=diff(XY(2:end));
    D=diff(XY);
    
    pause(0.01)
    subplot(2,2,2*fil-1)
    plot(XY(1:end-1,:),D, '.r');hold on;axis equal; 
    
    subplot(2,2,2*fil)
    plot(D1,D2, '.b');hold on;axis equal
end

end

%% Korrelation i TN
clc;clf;clear all

load('filnamn.mat')


for fil=1:2 %to plot data for both cell types

data =load(filnamn{fil});
C = separera(data);
n=length(C);%antal partiklar

figure(fil);clf


for i=1:n
    %if size(C{i},1)~=1000
    %    continue
    %end
    TN=koordinatbyte(C{i}(:,2:3));
    dTN=diff(TN(:),1,1);
    plot(dTN(1:end-1,1), dTN(2:end,1), '.b');hold on
end

end
















