%% filamentrörelser
clear all
D = dir('*.tif'); %hämtar alla filer i directory med formatet .tif
imcell = cell(1,numel(D));

%%
P = zeros(456,371,numel(D));
for i = 1:numel(D)
  imcell{i} = imread(D(i).name);
  P(:,:,i) = imcell{i};
end


%% "film" på strängen
bilder_s=20; %bilder/sekund
disp('speltid=')
disp(numel(D)/bilder_s)
for i = 1:numel(D)
  image(P(:,:,i))
  pause(bilder_s^-1)
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
A = importdata('confined_270304-6-28-min.mat');


