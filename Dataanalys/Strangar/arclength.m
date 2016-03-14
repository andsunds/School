function [ L ] = arclength( px, py, varargin )
%Beräknar båglängen på strängen.
%
%Om bara px och py ges, så antas längden vara från 0 till 1 och beräknas
%med 1000 punkters noggrannhet.
%
%Annars kan man ge start och slutpunkt också, man får då 1000 punkter
%noggrannhet.
%Till sista kan man även ge noggrannheten.


if nargin == 2
    S=linspace(0,1, 1000);
elseif nargin == 4
    S=linspace(varargin{1},varargin{2}, 1000);
elseif nargin == 5
    S=linspace(varargin{1},varargin{2}, varargin{3});
else
    fprintf('Fel antal inargument! \nFick %d, kan bara vara 2, 4 eller 5.', nargin);
    L=[];
    return
end


L=sum(sqrt(sum(diff([polyval(px, S); polyval(py, S)] ,2).^2 ,1)));

end

