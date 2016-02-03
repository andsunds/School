%% riktningsberoende

%% för ögoninspektion
clc;clf;clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=1;
data =load(filnamn{fil});

C = separera(data);
n=length(C);%antal partiklar


for i=1:n

    X=C{i}(:,2)-mean(C{i}(:,2));
    Y=C{i}(:,3)-mean(C{i}(:,3));
    M=max(max(abs([X,Y])));
    %Kollar fördelning av polära vinkelkoordinater
    %v=atan(Y./X);
    %v=atan(diff(Y)./diff(X));
    %figure(1)
    %hist(v)
    
    %Hittar min. kv. anpassning
    [koef, R_sq]=minsta_kvadrat( X, Y );
    
    figure(2)
    plot(X,Y, '.');hold on
    x=[-M, M]*1.1;
    lutn=-koef(1)/koef(2);
    y=lutn*x;
    p=plot(x,y);hold off
    axis([-M, M, -M, M]*1.1)
    axis equal
    legend(p,sprintf('y=%0.2fx, R^2=%0.2f',lutn,R_sq))
       
    pause(2)
end


%% Jämförelse mot storleken
clc;clf;clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';



for fil=1:2
    data =load(filnamn{fil});
    C = separera(data);
    n=length(C);%antal partiklar
    
    I=zeros(n,1);
    R_sq=zeros(n,1);
    
    for i=1:n
        I(i)=mean(C{i}(:,4));
        X=C{i}(:,2)-mean(C{i}(:,2));
        Y=C{i}(:,3)-mean(C{i}(:,3));
        [koef, R_sq(i)]=minsta_kvadrat( X, Y );
    end
    
    figure(fil)
    plot(I,R_sq,'.')
    title(filnamn{fil})
    xlabel('$I$', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
    ylabel('$R^2$', 'Interpreter', 'Latex', 'FontSize', 16, 'Color', 'k');
    set(gca,'FontSize',15)%,'XScale','log','YScale','log');
    figure(fil+2)
    hist(R_sq)
end








