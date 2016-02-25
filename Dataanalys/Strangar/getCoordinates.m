%% Plocka ut koordinater
clear all; clc;clf
filnamn=cell(1,4);
filnamn{1}='confined_270304-6-28-min.mat'; % Några konstiga hopp i denna.
filnamn{2}='confined_280204-2-32min.mat';
filnamn{3}='nonconfined_180304-1-5min.mat';
filnamn{4}='nonconfined_250104-1-167min.mat';

A = importdata(filnamn{4});
dim=2;
N_pics=size(A,3);

%Vi vill veta hur många punkter som finns för varje bild
N_points=zeros(N_pics,1);
for i = 1:N_pics
    N_points(i)=length(find(A(:,:,i)));
end

%Här hämtas koordinaterna

coordinates=NaN(N_pics, max(N_points), dim);

for i = 1:N_pics
    [row, kol]=find(A(:,:,i));
    %l=length(kol);
    coordinates(i,1:N_points(i), 1)=kol;
    coordinates(i,1:N_points(i), 2)=row;
end
%%
clf
L_string   = sum( sqrt( nansum( diff(coordinates,1, 2).^2 , 3 ) ), 2);

L_endtoend = zeros(N_pics,1);
for i = 1:N_pics
    L_endtoend(i) =... 
       sqrt( sum( ... 
         ( coordinates(i, N_points(i), :) - coordinates(i, 1, :) ).^2, ...
          3) );                         % ^avst. mellan första och sista punkten
        % ^summera längs tredje dimensionen (summera x och y)
end

subplot(1,2,1)
plot(L_string), hold on
plot(L_endtoend)
subplot(1,2,2)
plot(L_endtoend./L_string)

%% Ska BARA köras en gång!!!!
%%%save('nonconfined_167min.mat', 'coordinates', 'N_points', 'L_string', 'L_endtoend', '-mat')
%%%disp('save successfull')

%% 
clf
for i = 1:N_pics
plot(coordinates(i,:,1),coordinates(i,:,2), 'o-')
axis equal
pause(.1)
end

%%