%%
clc;clear all;clf

dim=3;               %dimension
N=3;m=ones(1,N); %m(1)=100000;m(end)=10000;

%m=[1000 1 1000];N=length(m); 

M=sum(m);            %mass of the particles

%k=[10; 100]; %spring constant
k=ones(N-1,1)*0; %spring constant
%R0=[.2; .5];%equilibrium distance
R0=ones(N-1,1); %equilibrium distance

s     = ones(N-2,1);
theta0= zeros(N-2,1);

damp=0;%dampening

%q0=[ -1,0,0,  0,.5,0, 1,0,0;     0,0,0,  0, 0,2, 0,0,0 ];

%No movement of CM
%q0=[ -.5 , 0 , -.3, .3, .1, -.3, 0, .7, 1,
%       -1 , 0 ,  0 , .7 , 1 , 0, 0, 0, 1, ];

%q0=[logspace(-1,1,N*dim);linspace(1,-1,N*dim)];
%q0(1,1)=0;        q0(1,2)=0;q0(2,1)=0;        q0(2,2)=0;
%q0(1,end-1)=-1.5; q0(1,end)=-3;q0(2,end-1)= 0  ; q0(2,end)=0;

%q0=randn([2,N*dim]);

q0=[ -1,0,0,  0,1,0, 1,0,0;     0,0,0,  0, 0,0, 0,0,0 ];
%q0=ones(2,N*dim).*[reshape(repmat(1:N,dim, 1),1,[]); zeros(1,N*dim)]

%q0=[1:N*dim; zeros(1,N*dim)];

play_movie_v2(0,reshape(q0,1,[]),m,dim,3,1 )




%q0=reshape(q0, 1,[])
q=reshape(q0, 2,[]);%makes a matrix from the vector

r=reshape(q(1,:), dim,[] ).';








dr=diff(r, 1, 1);

%radial force
k=reshape(k, [],1);
f_radial=dr.*k.*(1-R0./sqrt(sum(dr.^2, 2)));

N=size(r,1);
a=eye(N,N-1)-[zeros(1,N-1); eye(N-1)];


F_radial=a*f_radial;







%axial force
s=reshape(s, [],1);
mo=cross(dr(1:end-1,:),dr(2:end,:),2);

%sin_angle=sqrt(sum(mo.^2,2))./(sqrt(sum(dr(1:end-1,:).^2, 2)).*sqrt(sum(dr(2:end,:).^2, 2)));
%angle=asin(sin_angle)

cos_angle=dot(dr(1:end-1,:),dr(2:end,:),2)./(sqrt(sum(dr(1:end-1,:).^2, 2)).*sqrt(sum(dr(2:end,:).^2, 2)));
angle=acos(cos_angle);

normalize=@(A) A./sqrt(sum(A.^2, 2));

mo=normalize(mo);
Mo=mo.*angle.*s

f_axial=reshape([cross(dr(1:end-1,:),Mo,2),cross(dr(2:end,:),Mo,2)].', dim,[]).';
%f_axial=reshape(f_axial.', dim,[]).';

l=size(f_axial,1);

b=[repmat([1 -1 -1 1], N,1), zeros(N,l-1)];
for i=0:N-1
    b(1+i,:)=circshift(b(i+1,:), i*2,2);
end

b=b(:,4:end)

F_axial=b*f_axial;


%F_axial=b*f_axial;


F=a*f_radial+b*f_axial







%%

N=length(m);
a=eye(N,N-1)-[zeros(1,N-1); eye(N-1)];

F=a*f;

A=F./repmat(m, dim, 1).';

Dq=[q(2,:);
     reshape(A.', 1, [])-damp*q(2,:)];

Dq=reshape(Dq, [], 1);



















