%% Autokorrelation med Wiener-Khinchin teoremet
clf;clear all;clc;
filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

nbrlags = 999; % Välj 999 för steg och 1000 för position
for fil = 1:2
    AN = zeros(nbrlags,1);
    AT = zeros(nbrlags,1);
    data =load(filnamn{fil});
    C = separera(data);
    n=length(C);%antal partiklar
    np = 0; %Antalet partiklar som beräkningen använder
    
    
    for i=1:n
        if length(C{i})==1000 && C{i}(1,4) < 5
            np = np+1;
            TN=koordinatbyte(C{i}(:,2:3));
            T = TN(:,1);
            N = TN(:,2);
            dT = diff(T);
            dN = diff(N);


            Tw = fft(dT);
            Nw = fft(dN);

            autokorrT = ifft(abs(Tw).^2);
            autokorrN = ifft(abs(Nw).^2);
            autokorrT = autokorrT/autokorrT(1,1);
            autokorrN = autokorrN/autokorrN(1,1);

            AN = AN+autokorrT;
            AT = AT+autokorrN;
        end
    end
    AN = AN/np;
    AT = AT/np;
    
    
    if fil==1
        titel = 'Energydepleted';
    else
        titel = 'Logphasecells';
    end
    
    subplot(2,2,2*fil-1)
    plot(AT(1:end/2))
    legend('Tangentiell','Interpreter', 'Latex')
    title(titel,'Interpreter', 'Latex')
    subplot(2,2,2*fil)
    plot(AN(1:end/2))
    legend('Normal','Interpreter','Latex')
    title(titel,'Interpreter', 'Latex') 
end
set(0,'defaultAxesFontSize', 16)
