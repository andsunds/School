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

lutn_x=zeros(n,1); R_sq_v=zeros(n,1);
lutn_v=zeros(n,1); R_sq  =zeros(n,1);

for i=1:n

    X=C{i}(:,2)-mean(C{i}(:,2));
    Y=C{i}(:,3)-mean(C{i}(:,3));
    dX=diff(X);
    dY=diff(Y);
        
    %Hittar min. kv. anpassning
    [koef, R_sq(i)] =minsta_kvadrat( X, Y );
    [koef_v, R_sq_v(i)] =minsta_kvadrat( dX, dY );
    
    lutn_x(i)=-koef(1)/koef(2);
    lutn_v(i)=-koef_v(1)/koef_v(2);
    
%     subplot(1,2,1)
%     title('Rumskoordinat')
%     M=max(max(abs([X,Y])));
%     x=[-M, M]*1.1;
%     y=lutn_x(i)*x;
%     plot(X,Y, '.');hold on
%     p=plot(x,y);hold off; legend(p,sprintf('y=%0.2fx, R^2=%0.2f',lutn_x(i),R_sq(i)))
%     axis equal; axis([-M, M, -M, M]*1.1)
    
    
%     subplot(1,2,2)
%     title('Hastighet')
%     Mv=max(max(abs([dX,dY])));
%     vx=[-Mv, Mv]*1.1;
%     xy=lutn_v(i)*vx;
%     plot(dX,dY, '.');hold on
%     h=plot(vx,xy);hold off; legend(h,sprintf('y=%0.2fx, R^2=%0.2f',lutn_v(i),R_sq_v(i)))
%     axis equal; axis([-Mv, Mv, -Mv, Mv]*1.1)

    %pause(1)
end
disp('Färdig')


%% Korrelationmellan hastighet och position?
clc;clf
%subplot(1,2,1)
hist(lutn_x),%hist(lutn_v)
%plot(lutn_x,lutn_v, '.')
%subplot(1, 2, 2)
%hist(R_sq), hold on
%hist(R_sq_v)










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








