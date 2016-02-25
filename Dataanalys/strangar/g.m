function g = g(tau,s)
% Beräkning av Måns "g-funktion". Behöver troligtvis vektoriseras, väldigt 
%långsam för tillfället. 
% tau: Max skillnad i tid
%s: Max skillnad i (z1-z2)

g = zeros(tau,s+1);

prec = 100; % Noggrannhet


    for kk=0:tau
        k=kk/prec;
        for jj=-s/2:s/2
            j=jj/prec;
            f = @(p) exp(-(1+p.^2).*k+1i.*p.*j)./(1+p.^2);
            g(kk+1,jj+s/2+1) = abs(1/(2*pi).*integral(f,-Inf,Inf,'RelTol',0,'AbsTol',1e-2));
        end
    end
    
    T = linspace(0,tau,tau+1);
    S = linspace(-s/2,s/2,s+1);
    surf(S,T,g)

end

