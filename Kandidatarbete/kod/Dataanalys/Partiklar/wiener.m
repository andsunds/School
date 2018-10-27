%% Autokorrelation med Wiener-Khinchin teoremet
clf;clear all;clc;

load('filnamn.mat')


figure(1)
nbrlags = 999; % V�lj 999 f�r stegl�ngd och 1000 f�r position
for fil = 1:2
    AN = zeros(nbrlags,1);
    AT = zeros(nbrlags,1);
    data =load(filnamn{fil});
    C = separera(data);
    n=length(C);%antal partiklar
    np = 0; %Antalet partiklar som ber�kningen anv�nder
    
    
    for i=1:n
        if length(C{i})==1000 && C{i}(1,4) < 5 % Om man vill begr�nsa storlek
            np = np+1;
            TN=koordinatbyte(C{i}(:,2:3));
            T = TN(:,1); % Ta fram tangential och normalkoord
            N = TN(:,2);
            dT = diff(T); % Steg tangentiell
            dN = diff(N); % Steg normal


            Tw = fft(dT); 
            Nw = fft(dN);

            %Ta fram autokorrelation enligt Wiener-Khinchin
            %Normera med "perfekt korrelation"
            autokorrT = ifft(abs(Tw).^2); 
            autokorrN = ifft(abs(Nw).^2);
            autokorrT = autokorrT/autokorrT(1,1);
            autokorrN = autokorrN/autokorrN(1,1);
            
            % Summera ihop f�r alla partiklar
            AN = AN+autokorrT;
            AT = AT+autokorrN;
        end
    end
    %Normera med antalet partiklar s.a autkorro in [-1,1]
    AN = AN/np;
    AT = AT/np;
    
    
    if fil==1
        titel = 'Energydepleted';
    else
        titel = 'Logphasecells';
    end
    
    
    figure(1)
    subplot(2,2,2*fil-1)
    plot(AT(1:end/2))
    legend('Tangentiell','Interpreter', 'Latex')
    title(titel,'Interpreter', 'Latex')
    subplot(2,2,2*fil)
    plot(AN(1:end/2))
    legend('Normal','Interpreter','Latex')
    title(titel,'Interpreter', 'Latex') 
    figure(2)
    
    % Fouriertransformen av autokorrelationen
    ANtrans = fft(AN);
    ATtrans = fft(AT);
    
    Fs = 100;            % Sampling frequency
    T = 1/Fs;             % Sampling period
    L = 1000;             % Length of signal
    t = (0:L-1)*T;        % Time vector
    
    P2N = abs(ANtrans/L);
    P1N = P2N(1:L/2+1);
    P1N(2:end-1) = 2*P1N(2:end-1);
    
    P2T = abs(ATtrans/L);
    P1T = P2T(1:L/2+1);
    P1T(2:end-1) = 2*P1T(2:end-1);
    
    f = Fs*(0:(L/2))/L;
    
    subplot(2,2,2*fil-1)
    plot(f,P1T)
    legend('Tangentiell','Interpreter', 'Latex')
    title(titel,'Interpreter', 'Latex')
    subplot(2,2,2*fil)
    plot(f,P1N)
    legend('Normal','Interpreter','Latex')
    title(titel,'Interpreter', 'Latex') 
    
end
set(0,'defaultAxesFontSize', 16)
