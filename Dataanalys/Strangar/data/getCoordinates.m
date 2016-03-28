%% Plocka ut koordinater
clear all; clc;clf

lengthscale=64.5e-9; % m/px, storlek av en pixel i bilden i meter


filnamn=cell(1,4);
filnamn{1}='confined_270304-6-28-min.mat'; % Några konstiga hopp i denna.
filnamn{2}='confined_280204-2-32min.mat';
filnamn{3}='nonconfined_180304-1-5min.mat';
filnamn{4}='nonconfined_250104-1-167min.mat';

fil=4;

a = importdata(filnamn{fil});

A=zeros((size(a)+[4, 4, 0]));
%disp('hej')
A(3:end-2, 3:end-2,:)=a;

dim=2;
N_pics=size(A,3);

%Vi vill veta hur många punkter som finns för varje bild
N_points=zeros(N_pics,1);
for i = 1:N_pics
    N_points(i)=length(find(A(:,:,i)));
end

%Här hämtas koordinaterna

coordinates=NaN(N_pics, max(N_points), dim);
index_endpoint=zeros(N_pics,1);

for i = 1:N_pics
    [rad, kol]=find(A(:,:,i));
    index_endpoint(i)=find_endpoint(A(:,:,i) , rad, kol );
    %l=length(kol);
    coordinates(i,1:N_points(i), 1)=kol;
    coordinates(i,1:N_points(i), 2)=rad;
    
    
end


coordinates=coordinates*lengthscale;


%% Sort points along the string
clf;clc
i=50;
S=coordinates;

%143
for i=1:N_pics %Över alla bilder
    %Lista med lediga index
    INDEX=1:N_points(i); 
    INDEX(index_endpoint(i))=[];
    %Initierar Smed första 
    S(i,1,:)=coordinates(i,index_endpoint(i),:);
    
    for j=2:N_points(i) %Över alla punkter
        [S(i,j,:), I] =find_nearest_point(S(i,j-1,:), coordinates(i,INDEX,:));
        %I
        INDEX(I)=[];
        
    end
end


coordinates=S;%sorterat!

% %% film
% clf
% for i = 1:N_pics
%     plot(S(i,:,1),S(i,:,2), '-'), hold on
%     plot(S(i,1,1),S(i,1,2), 'r*', 'markersize', 16), hold off
%     pause(.1)
% 
% end

%% Beräkna längd
clf;clc

L=cumsum( sqrt( nansum( diff(coordinates, 1, 2).^2 , 3 ) ), 2);

L_string=[zeros(size(L,1),1), L];

L_endtoend = zeros(N_pics,1);
for i = 1:N_pics
    L_endtoend(i) =... 
       sqrt( sum( ... 
         ( coordinates(i, N_points(i), :) - coordinates(i, 1, :) ).^2, ...
          3) );                         % ^avst. mellan första och sista punkten
        % ^summera längs tredje dimensionen (summera x och y)
end
size(L_string)
size(L_endtoend)

subplot(1,2,1)
plot(L_string(:,end)), hold on
plot(L_endtoend)
%axis([0 N_pics min(L_endtoend)-50 max(max(L_string))+50])
title('Length and end-to-end length')
xlabel('pictureframe')
legend('Length','End-to-end','Location','Best')

subplot(1,2,2)
plot(L_endtoend./L_string(:,end))
%axis([0 N_pics min(L_endtoend./L_string(:,end))-0.05 1])
title('end-to-end/length')
xlabel('pictureframe')


% %%
% clc;clf
% R_sq=mean((L_endtoend.^2));
% L=mean(L_string);
% 
% f=@(x) R_sq -2*(L.*x + x.^2.*(-1 + exp(-L./x)));
% 
% x=linspace(0,1e3);
% plot(x,f(x))
% 
% fzero(f, 100)

%% Ska BARA köras en gång!!!!
filnamn=cell(1,4);
filnamn{1}='confined_28min.mat'; % Några konstiga hopp i denna.
filnamn{2}='confined_32min.mat';
filnamn{3}='nonconfined_5min.mat';
filnamn{4}='nonconfined_167min.mat';


save(filnamn{fil}, 'coordinates', 'N_points', 'L_string', 'L_endtoend', '-mat')
disp('save successfull')
