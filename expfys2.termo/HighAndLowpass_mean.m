function [xl,xh] = HighAndLowpass_mean(t,x,T)
%HIGHPASS_MEAN Högpassfiltrerar med ett löpande medelvärde över en period.

t = t(:,1);
xl = zeros(size(x));
xh = zeros(size(x));
for i=1:size(x,2)
    % Första och sista perioden används bara medelvärdet för den
    % första/sista perioden. Annars får man systematisk avvikelse.
    
    % Observera att first:last och 1:f överlappar
    first = find(t-t(1) >= T/2, 1, 'first');
    last = find(t(end)-t >= T/2, 1, 'last');
    f = find(t-t(1) >= T, 1, 'first')-1;
    l = find(t(end)-t >= T, 1, 'last')+1;
    
    tmp = smooth(t,x(:,i),T,'moving');
    xl(1:f,i) = trapz(t(1:f),x(1:f,i))/(t(f)-t(1)); % Medelvärde av början
    xl(l:end,i) = trapz(t(l:end),x(l:end,i))/(t(end)-t(l)); % Medelvärde av slutet
    xl(first:last,i) = tmp(first:last); % Lågpass av mitten
    
    xh(:,i) = x(:,i)-xl(:,i); %Högpass = signal-lågpass
end

end
