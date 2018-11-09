%%
clc;clf;clear

data_verlet=load('verlet.tsv');

t_v=data_verlet(:,1);
x_v=data_verlet(:,2);
v_v=data_verlet(:,3);
E_v=0.5*v_v.^2+0.5*x_v.^2;

plot(t_v,E_v)
%plot(t_v,x_v)


%%
clc;clf;clear
data_euler=load('euler2F.tsv');

t_e=data_euler(:,1);
x_e=data_euler(:,2);
v_e=data_euler(:,3);
E_e=0.5*v_e.^2+0.5*x_e.^2;

hold on
plot(t_e,E_e)


dt=mean(diff(t_e));
w=1;

%T=E_e(1)*exp(dt*t_e/w^2);
%plot(t_e,T)




