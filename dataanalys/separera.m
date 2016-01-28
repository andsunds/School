function [C,index] = separera(X) % Separera partiklar och spara i en cell

% B = [1;X(:,4)];
% P = [X(:,4);1];
% E = B-P;
% index = find(E);


index=[find(X(:,1)==0); size(X,1)+1];%varje index som en ny partikel börjar vid

n = size(index,1)-1; % n �r antalet olika partiklar

C = cell(n,1);
for i=1:n
    C{i} = X(index(i):index(i+1)-1, :);
end



