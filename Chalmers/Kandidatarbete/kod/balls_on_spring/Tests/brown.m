%% Single particle Brwonian motion (cumsum)
clc;clf;clear all

N=1000; 

tic
Hx=rand(N,1)-.5;
Hy=rand(N,1)-.5;

x=cumsum(Hx);
y=cumsum(Hy);
toc

% fast plot
subplot(1,2,1)
plot(x, y, '.-');

subplot(1,2,2)
d=sum([x,y].^2, 2);
j=0:(N-1);
plot(j,d);



% slow plot

%distance
subplot(1,2,2)
d=sqrt(sum([x,y].^2, 2));
j=0;
D=d(1);
p=plot(j,D);
p.XDataSource = 'j';
p.YDataSource = 'D';



%position
subplot(1,2,1)
X=x(1);Y=y(1);
h=plot(X, Y, '.-');
grid on
hold on

M=max(abs([x;y]));

h.XDataSource = 'X';
h.YDataSource = 'Y';


for i=1:N
    X=x(1:i);Y=y(1:i);
    j=0:(i-1); D=d(1:i);
    refreshdata
    axis([-M, M, -M, M]*1.1)
    pause(0.001)
end


%% Single particle Brwonian motion (ODE45)
clc;clf;clear all
%this give increasing solutions
r0=[0,0];
odefun=@(x, t) rand(2,1);
[T, R] = ode45(odefun, [0,100], r0);

%Compensate with '-.5*T'
x=R(:,1)-.5*T;
y=R(:,2)-.5*T;
%plot(x,y)


N=length(T);

X=x(1);Y=y(1);
h=plot(X, Y, '.-');
grid on
hold on
axis equal

M=max(abs([x;y]));

h.XDataSource = 'X';
h.YDataSource = 'Y';

dt=diff(T);

for i=1:N
    X=x(1:i);Y=y(1:i);
    refreshdata
    %axis([-M, M, -M, M]*1.1)
    pause(dt(i))
end









































