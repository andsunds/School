function [ koef , R_sq ] = minsta_kvadrat( X, Y )
%Minsta kvadratanpassning där det vinkelräta avståndet till 
%linjen minimeras.
%
%Viktigt att X och Y vardera har medelvärde 0.


[koef, SS_res]=fminsearch(@(L) f(L(1),L(2), X,Y), [1;1]);

SS_tot=sum(X.^2+Y.^2);

R_sq=1-SS_res/SS_tot;

end


function [dist_sq] = f(a, b, x,y)
%Beräknar summan av minsta vinkelräta avståndet i kvadrat 
%från linjen.

dist_sq = sum((a*x +b*y).^2./(a.^2+b.^2));
end
