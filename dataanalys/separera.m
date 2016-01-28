function [C,index] = separera(X) % Separera partiklar och spara i en cell



index=[find(X(:,1)==0); size(X,1)+1];%varje index som en ny partikel börjar vid

n = size(index,1)-1; % n = antalet olika partiklar

C = cell(n,1);
for i=1:n
    C{i} = X(index(i):index(i+1)-1, :);
end



