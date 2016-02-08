%% riktningsberoende

%% för ögoninspektion utan koordinatbyte
clc;clf;clear all

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=2;
data =load(filnamn{fil});

C = separera(data);
n=length(C);%antal partiklar


for i=1:n

    X=C{i}(:,2)-mean(C{i}(:,2));
    Y=C{i}(:,3)-mean(C{i}(:,3));
    dX=diff(X);
    dY=diff(Y);
        
    %Hittar min. kv. anpassning
    [koef, R_sq] =minsta_kvadrat( X, Y );
    [koef_d, R_sq_d] =minsta_kvadrat( dX, dY );
    
    figure(1)
    subplot(1,2,1)
    title('Rumskoordinat')
    plot(X,Y, '.');hold on
    M=max(max(abs([X,Y])));
    x=[-M, M]*1.1;
    lutn=-koef(1)/koef(2);
    y=lutn*x;
    p=plot(x,y);hold off
    axis equal
    axis([-M, M, -M, M]*1.1)
    legend(p,sprintf('y=%0.2fx, R^2=%0.2f',lutn,R_sq))
    
    subplot(1,2,2)
    title('Hastighet')
    plot(dX,dY, '.');hold on
    M=max(max(abs([dX,dY])));
    x=[-M, M]*1.1;
    lutn=-koef_d(1)/koef_d(2);
    y=lutn*x;
    h=plot(x,y);hold off
    axis equal
    axis([-M, M, -M, M]*1.1)
    legend(h,sprintf('y=%0.2fx, R^2=%0.2f',lutn,R_sq_d))
    
    
    
    pause(.1)
end



%% för ögoninspektion med koordinatbyte
clc;clf;clear all
grid on

filnamn=cell(1,2);
filnamn{1}='energydepletedcells.csv';
filnamn{2}='logphasecells.csv';

fil=1;
data =load(filnamn{fil});

C = separera(data);
n=length(C);%antal partiklar


for i=1:n

    TN=koordinatbyte( bsxfun(@minus, C{i}(:,2:3), mean(C{i}(:,2:3), 1)) );
    
    T=TN(:,1);
    N=TN(:,2);
    l=floor(length(T)/2);
    
    figure(3)
    T_trans=fft(diff(T));
    loglog(abs(T_trans(2:l)));grid on
    figure(4)
    N_trans=fft(diff(N));
    loglog(abs(N_trans(2:l)));grid on
    
    
    
    figure(2)
    plot(T,N, '.');grid on,%hold on
    M=max(max(abs([T,N])));
    %t=[-M, M]*1.1;
    %lutn=-koef(1)/koef(2);
    %n=lutn*t;[n(2); -n(1)];

    %p=plot(t,n);hold off
    axis equal
    axis([-M, M, -M, M]*1.1)
    
    %legend(p,sprintf('y=%0.2fx, R^2=%0.2f',lutn,R_sq))
       
    pause(.1)
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
        X=diff(C{i}(:,2)-mean(C{i}(:,2)));
        Y=diff(C{i}(:,3)-mean(C{i}(:,3)));
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








