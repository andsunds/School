function [ koef ] = storleksanpassning( fil )
%Tar fram potenssamband för att vikta om ett medelvärde beroende på
%partikelns intensitet.

load('../filnamn.mat')

load(['../', kompl],...
      'intensitet', 'std_n', 'std_t', 'sigma_brus', '-mat')

I=intensitet{fil};

%Detta är just nu samma normering som för rörligheten...
lambda=sqrt(std_n{fil}.^2+std_t{fil}.^2 -2*sigma_brus(fil).^2 );
%beräkna koefficienter för lutning i loglog
koef=[ones(size(I)), log(I)]\log(lambda);
koef(1)=exp(koef(1));%konverterar till potenssamband

end

