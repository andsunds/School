%% filamentr�relser
clear all
D = dir('*.tif'); %h�mtar alla filer i directory med formatet .tif
imcell = cell(1,numel(D));

%%
%P = zeros(456,371,numel(D));
for i = 1:numel(D)
  imcell{i} = imread(D(i).name);
  %P(:,:,i) = imcell{i};
end

%% trans. hast.
t=5*60/numel(D);    %tid om hela f�rloppet �r 5minuter
L=108*10^-9;        %om 1pixel motsvarar 108nm
for i=1:numel(D)
[rad,kol]=find(255==imcell{i});
R(i)=sum(rad)./length(rad);
K(i)=sum(kol)./length(kol);
if i>1
    hast(i)=L/t*norm([R(i)-R(i-1),K(i)-K(i-1)]); %filementets hastighet
end
end

%% 
clear all
filnamn=cell(1,4);
filnamn{1}='confined_270304-6-28-min.mat'; % N�gra konstiga hopp i denna.
filnamn{2}='confined_280204-2-32min.mat';
filnamn{3}='nonconfined_180304-1-5min.mat';
filnamn{4}='nonconfined_250104-1-167min.mat';

A = importdata(filnamn{4});
% "film" p� str�ngen
bilder_s=20; %bilder/sekund
disp('speltid=')
disp(size(A,3)/bilder_s)
for i = 1:size(A,3)
  image(A(:,:,i))
  pause(bilder_s^-1)
end
%%
clc, clear all
A=importdata('nonconfined_250104-1-167min.mat');
for t=1:size(A,3)
[rad,kol]=find(255==A(:,:,t));
AA=A(:,:,t);
for i=1:length(rad)
    for j=1:5
        for k=1:5
    B(k,j)=A(rad(i)+k-3,kol(i)+j-3,t);
    [radB,kolB]=find(255==B);
    C(i)=length(radB);
        end
    end
end
D=find(3==C);
reserv=find(4==C);
if length(D)==1
    D=reserv;
end
if norm(D)==0
    D=reserv;
end

startXY(t,:)=[kol(D(1)),rad(D(1))];
%tail=[rad(D(2)),kol(D(2))];

%plot(kol,rad,start(2),start(1),'o')
%pause(.3)
end



