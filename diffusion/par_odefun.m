function [ Df ] = par_odefun(f, tau, alpha, dx)

%df=zeros(size(f));

d2f=diff(f,2);

df0=d2f(2)-d2f(1);
dfend=d2f(end-1)-d2f(end);


Df=alpha*[df0*0; d2f ; 0*dfend]/dx;




end
