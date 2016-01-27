%% model of two balls on a spring
clc;clf;clear;

m1=inf; m2=1; M=m1+m2; %mass of the two particles
k=10;                 %spring constant
n=1;                 %dampening
R0=0.25;              %equilibrium distance
Br=0;                %Browninan impact strength

v0x=.1; %initial speed
v0y=.5; 

%No movement of CM
q0=[ -.5 , 0 , .5, 0
    -v0x/m1, v0y/m1, v0x/m2, -v0y/m2];
  
odefun=@(t, q) dq_two_balls(q, m1, m2, k, n, R0, Br);
%reshape(q0, [],8) %test
%odefun(0, q0)'%test
 
[T, Q] = ode45(odefun, [0,100], reshape(q0,[], 8));


%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%
dt=diff(T);
N=length(T);

%position of particles
x=[Q(:,1) Q(:,5)];
y=[Q(:,3) Q(:,7)];

%CM pos.
x_cm=(m1*x(:,1)+m2*x(:,2))/M;
y_cm=(m1*y(:,1)+m2*y(:,2))/M;

%Individual particle
x1=x(:,1); y1=y(:,1);
X1=x1(1) ; Y1=y1(1);
x2=x(:,2); y2=y(:,2);
X2=x2(1) ; Y2=y2(1);



p1=plot(X1,Y1,'ob'); 
hold on
axis equal
p1.XDataSource = 'X1';
p1.YDataSource = 'Y1';

p2=plot(X2,Y2,'or'); 
p2.XDataSource = 'X2';
p2.YDataSource = 'Y2';


X=x(1,:);
Y=y(1,:);
p3=plot(X,Y,'-k');
p3.XDataSource = 'X';
p3.YDataSource = 'Y';


X_cm=x_cm(1);
Y_cm=y_cm(1);
p4=plot(X_cm,Y_cm,'xr');
p4.XDataSource = 'X_cm';
p4.YDataSource = 'Y_cm';



%movietime!
for i=1:N
    X1=x1(i);     Y1=y1(i);
    X2=x2(i);     Y2=y2(i);
    X=x(i,:);     Y=y(i,:);
    X_cm=x_cm(i); Y_cm=y_cm(i);
    refreshdata
    axis([-1 1 -1 1])
    pause(dt(i))
end
















