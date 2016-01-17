function [  ] = test_movie( f, t_start, t_end )

t=linspace(t_start,t_end);
dt=[0,diff(t)];
%global T Y p
T=t(1);Y=f(T);
p=plot(T,Y);hold on
h=plot(T,f(T/2));

%p.XDataSource = 'T';
%p.YDataSource = 'Y';

for i=1:length(t)
    pause(dt(i))
    T=t(1:i);Y=f(T);
    %global T Y
    %refreshdata(p)
    set(p,'XData',T,'YData',Y);set(h,'XData',T,'YData',f(T/2));
end

