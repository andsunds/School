function xf = highpass(t,x,fc)
%HIGPASS Funktion som med fft högpassfiltrerar en signal
% If-satserna är för att det bara funkar med vektorer av udda längd

N = length(x);
M = floor((N-1)/2);

if 2*M+1 ~= N
    % N är jämn
    x0 = x(2:end,:);
    N = N-1;
else
    x0 = x;
end

y = fft(x0);
y_tilde = fftshift(y);

w = -M:M;
w = w*2*pi/(t(end)-t(1));
wc = 2*pi*fc;
F = 1./(1+1i*wc./w);

yf_tilde = bsxfun(@times,y_tilde,reshape(F,N,1));
yf = ifftshift(yf_tilde);
xf = ifft(yf);

if 2*M+1 ~= length(x)
    xf = [x(1,:) ; xf];
end

end
