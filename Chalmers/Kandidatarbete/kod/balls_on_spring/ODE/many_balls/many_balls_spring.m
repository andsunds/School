%% model of N balls on a spring
clc;clf;clear;

dim=2;               %dimention
N=20;

%m=[1 4 .5, 2, 1];
m=ones(1,N);
m(1)=100000;m(end)=100;

M=sum(m);            %mass of the particles
%N=length(m);         %number of particles
%k=[10; 100; 20; 2];  %spring constant
k=ones(N-1,1)*100;  %spring constant
n=0;                 %dampening
%R0=[.2; .5; .1; 0];              %equilibrium distance
R0=ones(N-1,1)*.1;              %equilibrium distance
Br=0;                %Browninan impact strength



%No movement of CM
%q0=[ -.5 , 0 , -.3, .3, .1, -.3, 0, .7, 1, 1
%       -1 , 0 ,  0 , .7 , 1 , 0, 0, 0, 1, -.2  ];
q0=[logspace(-1,0,N*dim);linspace(1,-1,N*dim)];

%q0(1,floor(N*7/8))=1;
%q0(1,floor(N/2)+1)=1;
%q0(1,floor(N/2)+2)=2;
q0(1,1)=0;        q0(1,2)=0;
q0(2,1)=0;        q0(2,2)=0;
q0(1,end-1)=-1.5; q0(1,end)=-3;
q0(2,end-1)= 0  ; q0(2,end)=0;






odefun=@(t, q) dq_many(q, m, k, n, R0, Br, N, dim);
%reshape(q0, [],8) %test
%odefun(0, q0)'%test

l=sum(size(q0));


[T, Q] = ode45(odefun, [0,100], reshape(q0,1, []));


%%%%%%%%%%%%%% plotting %%%%%%%%%%%%%%
dt=diff(T);
L=length(T);


%position of particles
step=2*dim;
x=Q(:, 1:step:end);
y=Q(:, 3:step:end);
%size(x)
m_matrix=repmat(m, L,1);

%CM pos.
x_cm=sum(m_matrix.*x, 2)/M;
y_cm=sum(m_matrix.*y, 2)/M;

%plots the bonds
X=x(1,:);
Y=y(1,:);
p=plot(X,Y,'-k'); hold on
p.XDataSource = 'X';
p.YDataSource = 'Y';

%different colors for each particle
%c=logspace(0,-3, N).';
c=linspace(1, 0, N).';
C=[c, circshift(c, [1 0]), circshift(c, [2 0])]*0;
s=scatter(X,Y, 81, C);
s.XDataSource = 'X';
s.YDataSource = 'Y';

%shows center of mass
X_cm=x_cm(1);
Y_cm=y_cm(1);
h=plot(X_cm,Y_cm,'xr');
h.XDataSource = 'X_cm';
h.YDataSource = 'Y_cm';



axis equal
grid on


% movietime!
R=norm(R0)*12; %frame size
%R=2;
for i=1:L
    %i
    %X1=x1(i);     Y1=y1(i);
    %X2=x2(i);     Y2=y2(i);
    X=x(i,:);     Y=y(i,:);
    X_cm=x_cm(i); Y_cm=y_cm(i);
    refreshdata
    axis([X_cm-R, X_cm+R, Y_cm-R, Y_cm+R])
    %pause(0.001)
    pause(dt(i))
end





































