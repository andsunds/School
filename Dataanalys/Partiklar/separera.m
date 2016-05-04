function [C,index] = separera(X) % Separera data fran olika partiklar och spara i en cell

X(:,2:3) = X(:,2:3)*1e-6; %Korrigera enhet l�ngder
%X(:,1) = X(:,1); %Korrigera tidsenhet 

index=[find(X(:,1)==0); size(X,1)+1];%Tar ut varje index som en ny partikel startar vid, letar efter t==0, tar med slutindex+1 pga for-loopen nedan

n = size(index,1)-1; % n = antalet olika partiklar, drar bort 1 pga tillagget ovan

C = cell(n,1); %Spara enskilda partiklars data i cellen
for i=1:n
    C{i} = X(index(i):index(i+1)-1, :); %-1 eftersom andra partikelns initiala datapunkt inte ska komma med
end



