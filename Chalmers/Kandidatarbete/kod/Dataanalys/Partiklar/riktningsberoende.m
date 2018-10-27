%% riktningsberoende

%% för ögoninspektion utan koordinatbyte
clc;clf;clear all
load('filnamn.mat')

fil=2;
data =load(filnamn{fil});

C = separera(data);
n=length(C);%antal partiklar

lutn1=zeros(n,1); R_sq=zeros(n,1);
kvot=zeros(n,1);
egen=zeros(n,2);

for i=1:n

    X=C{i}(:,2)-mean(C{i}(:,2));
    Y=C{i}(:,3)-mean(C{i}(:,3));
    
    
    gyr = [sum(X.^2), sum(X.*Y); 
           sum(Y.*X), sum(Y.^2)]/length(X);
    
    %N_steps=length(X);
    %gyr = [sum(X.^2)-sum(X).^2/N_steps, sum(X.*Y)-sum(X).*sum(Y)/N_steps; 
    %       sum(Y.*X)-sum(X).*sum(Y)/N_steps, sum(Y.^2)-sum(Y).^2/N_steps]/N_steps;
    [V,D]=eig(gyr);
    kvot(i)=D(4)/D(1);
    egen(i,:)=[D(1), D(4)];
        
    %Hittar min. kv. anpassning
    [koef, R_sq(i)] =minsta_kvadrat( X, Y );
    
    
    lutn1(i)=-koef(1)/koef(2);
    lutn2=V(2,2)/V(1,2);
    
    %{
    title('Rumskoordinat')
    M=max(max(abs([X,Y])));
    x=[-M, M]*1.1;
    y1=lutn1(i)*x;
    y2=lutn2*x;
    plot(X,Y, '.');hold on
    p=plot(x,y1);
    q=plot(x,y2);
    hold off; 
    legend([p,q],...
           sprintf('y=%0.2fx, R^2=%0.2f',lutn1(i),R_sq(i)),...
           sprintf('y=%0.2fx',lutn2) )
    axis equal; axis([-M, M, -M, M]*1.1)
    
    
    pause(1)
    %}
end
bins=50;
A=(egen(:,2)-egen(:,1)).^2./(egen(:,2)+egen(:,1)).^2;%assymetri
    [N_X, X]=hist(A,bins);
    hold on
    plot(X+.5*[diff(X),0], 1-cumsum(N_X)/sum(N_X), '-o')


disp('Färdig')

% % save('simuleringar/kvot_logphase.mat', 'kvot', 'egen', '-mat')


%% Korrelationmellan hastighet och position?
clc;clf
%subplot(1,2,1)
hist(lutn1),%hist(lutn_v)
%plot(lutn_x,lutn_v, '.')
%subplot(1, 2, 2)
%hist(R_sq), hold on
%hist(R_sq_v)






%% Jämförelse mot storleken
clc;clf;clear all
load('filnamn.mat')



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








