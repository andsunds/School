%% filamentrörelser
clear all
D = dir('*.tif'); %hämtar alla filer i directory med formatet .tif
imcell = cell(1,numel(D));

%%
%P = zeros(456,371,numel(D));
for i = 1:numel(D)
  imcell{i} = imread(D(i).name);
  %P(:,:,i) = imcell{i};
end

%% trans. hast.
t=5*60/numel(D);    %tid om hela förloppet är 5minuter
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
filnamn{1}='confined_270304-6-28-min.mat'; % Några konstiga hopp i denna.
filnamn{2}='confined_280204-2-32min.mat';
filnamn{3}='nonconfined_180304-1-5min.mat';
filnamn{4}='nonconfined_250104-1-167min.mat';

A = importdata(filnamn{1});
% "film" på strängen
bilder_s=20; %bilder/sekund
disp('speltid=')
disp(size(A,3)/bilder_s)
for i = 1:size(A,3)
  image(A(:,:,i))
  pause(bilder_s^-1)
end
