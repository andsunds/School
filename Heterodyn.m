function [ log_r, phi ] = Heterodyn(time,volt,period,nper,order)

if size(time,2) == 2
    time = time(:,1) + 0.5*diff(time,1,2);
    time = repmat(time,1,size(volt,2));
elseif size(time,2)>size(volt,2)
    time = time(:,1:end-1)+.5*diff(time,1,2);
else
    t=reshape(time.', [], 1);
    dt=diff(t);
    time= reshape(t+.5*[dt; mean(dt)], [],size(volt,1)).';
end

Z = zeros(size(time,2),1);
for k=1:size(time,2);
    n = find(time(end,k)-time(:,k) < nper*period,1);

    time_k = time(n:end,order(k));
    volt_k = volt(n:end,order(k));
    % Interpolera med styckvis 3:e-gradspolynom
    pp = spline(time_k,volt_k); %pp=piecewise polynomial
    % Sampla interpolationen oftare
    t = (time_k(1):(period*1e-3):time_k(end)).';
    % Multiplicera med kernel
    P = exp(-1i*2*pi/period*t).*ppval(pp, t);
    % Integrera
    I = trapz(t,P,1);
    % Medelv�rdera (se f�rstudie f�r 2*xxx)
    Z(k) = 2*I./(t(end)-t(1));
end
%I = sum(P(1:end-1,:).*diff(time,1,1));
%disp(abs(Z))

log_r = log(abs(Z));
phi = phase(Z);
for i=2:length(phi)
    if(phi(i) > phi(i-1))
        phi(i) = phi(i) - 2*pi*ceil((phi(i)-phi(i-1))/(2*pi));
    end
end

end
