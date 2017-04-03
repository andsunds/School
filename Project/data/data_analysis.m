%% Some plots 
% Plots E and M as a function of time
clc;clear all;clf

FILENAME = 'Ising/EM_beta_1.00000_PERIODIC.bin';
fID=fopen(FILENAME,'rb');
data1=fread(fID,[2,inf],'real*8').';
fclose(fID);
FILENAME = 'Ising/EM_beta_3.00000_PERIODIC.bin';
fID=fopen(FILENAME,'rb');
data3=fread(fID,[2,inf],'real*8').';
fclose(fID);

figure(5)

subplot(1,2,1)
plot(data1(:,1), 'b-'), hold on
plot(data3(:,1),'r-.')
set(gca, 'xlim',[0,1e7], 'xscale','log', 'fontsize',13)
grid on
xlabel('$i$', 'interpreter', 'LaTeX')
ylabel('$E_i$', 'interpreter', 'LaTeX')
l=legend('$J\beta=1$','$J\beta=3$');
set(l, 'interpreter', 'LaTeX')

subplot(1,2,2)
plot(data1(:,2), 'b-'), hold on
plot(data3(:,2),'r-.')
set(gca, 'yLim',[-1,1], 'xlim',[0,1e7], 'xscale','log', 'fontsize',13)
grid on
xlabel('$i$', 'interpreter', 'LaTeX')
ylabel('$M_i$', 'interpreter', 'LaTeX')
l=legend('$J\beta=1$','$J\beta=3$');
set(l, 'interpreter', 'LaTeX')

clear data1 data3



%% Ising
clc;clear all

%need to redo a sweep through the spectrum


data = load('Ising/EstdEMstdM_beta_0.100-2.147_2048_PERIODIC.tsv', '-ascii');



N=16*16;

T = data(:,1);
E = data(:,2)/N; stdE = data(:,3);%/sqrt(2e7);
M = data(:,4);   stdM = data(:,5);%/sqrt(2e7);

N_mean=32;
TT=   mean(reshape(T, N_mean,[]),1).';
EE=   mean(reshape(E, N_mean,[]),1).';
stdEE=mean(reshape(stdE, N_mean,[]),1).';
MM=   mean(reshape(abs(M), N_mean,[]),1).';
stdMM=mean(reshape(stdM, N_mean,[]),1).';


figure(1), clf
subplot(1,2,1)
plot(TT,EE, 'b.')


set(gca, 'xScale', 'lin', 'fontsize', 13)
set(gca, 'ylim', [-2, 0]);
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$\langle E\rangle/(JN)$', 'interpreter', 'LaTeX');


subplot(1,2,2)
plot(TT,abs(MM), 'b.')
%plot(TT,stdMM, '.')

set(gca, 'yScale', 'lin', 'xscale','lin', 'fontsize', 13)
set(gca, 'xlim',[0,4])
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$|\langle M\rangle|$', 'interpreter', 'LaTeX');

%% 1C: CV
clc,clear all

data2 = load('Ising/EstdEMstdM_beta_0.350-0.555_8192_PERIODIC.tsv');
data = load('Ising/EstdEMstdM_beta_0.100-2.147_2048_PERIODIC.tsv', '-ascii');

N=16*16;

T = data(:,1); T2 = data2(:,1); 
stdE = data(:,3);
stdE2 = data2(:,3);


N_mean=16;
N_mean2=128;
TT=   mean(reshape(T, N_mean,[]),1).';
TT2=   mean(reshape(T2, N_mean2,[]),1).';
stdEE=mean(reshape(stdE, N_mean,[]),1).';
stdEE2=mean(reshape(stdE2, N_mean2,[]),1).';


figure(2), clf
CV=stdEE.^2./TT.^2;
CV2=stdEE2.^2./TT2.^2;

hold on
plot(TT2,CV2/max(CV2),'b.')
plot(TT(1:64),CV(1:64)/max(CV(21)),'rx')

set(gca, 'xScale', 'lin', 'fontsize', 13)
set(gca, 'xlim', [0, 6]);
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$\hat{C}_V$', 'interpreter', 'LaTeX');




%% 1 D
clc;clear all

%need to redo a sweep through the spectrum

data12 = load('Ising/EstdEMstdM_T_1.400-3.447_2048_L-12_PERIODIC.tsv');
data16 = load('Ising/EstdEMstdM_T_1.400-3.447_2048_L-16_PERIODIC.tsv');
data20 = load('Ising/EstdEMstdM_T_1.400-3.447_2048_L-20_PERIODIC.tsv');
data24 = load('Ising/EstdEMstdM_T_1.400-3.447_2048_L-24_PERIODIC.tsv');

L=12:4:24;

Nmean=8;
TT=data12(:,1);
T = mean(reshape(TT, Nmean,[]), 1).';
chi12=L(1)^2*mean(reshape(data12(:,5).^2./TT, Nmean,[]),1).';
chi16=L(2)^2*mean(reshape(data16(:,5).^2./TT, Nmean,[]),1).';
chi20=L(3)^2*mean(reshape(data20(:,5).^2./TT, Nmean,[]),1).';
chi24=L(4)^2*mean(reshape(data24(:,5).^2./TT, Nmean,[]),1).';




figure(3), clf

hold on
plot(T,chi12, 'b*')
plot(T,chi16, 'rx')
plot(T,chi20, 'k.')
plot(T,chi24, 'm+')

set(gca, 'xScale', 'lin', 'fontsize', 13)
%set(gca, 'ylim', [-2, 0]);
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$\chi$', 'interpreter', 'LaTeX');
l=legend('$L=12$','$L=16$','$L=20$','$L=24$');
set(l, 'interpreter', 'LaTeX')

figure(4), clf
hold on
Tc=2/log(1+sqrt(2));
t=(T-Tc)/Tc;

nu=1.1; a=1/nu;
gamma=1.7; b=gamma/nu;

%a = 1.2;
%b = 1.75; 
fprintf('nu = %1.2f\n', 1/a);
fprintf('gamma = %1.2f\n', b/a);



tt  = repmat(t,1,4).*repmat(L.^(a),length(t),1);
chi = [chi12, chi16, chi20, chi24].*repmat(L.^(-b),length(t),1);

plot(tt(:,1),chi(:,1), 'b--')
plot(tt(:,2),chi(:,2), 'r:')
plot(tt(:,3),chi(:,3), 'k-')
plot(tt(:,4),chi(:,4), 'm-.')

%plot([0,0],[0,max(max(chi))],'-k', 'linewidth',2)

axis([-5,10, 0,1])

grid on
set(gca, 'xScale', 'lin', 'fontsize', 13)
xlabel('$tL^{1/\nu}$', 'interpreter', 'LaTeX');
ylabel('$\chi L^{-\gamma/\nu}$', 'interpreter', 'LaTeX');
l=legend('$L=12$','$L=16$','$L=20$','$L=24$');
set(l, 'interpreter', 'LaTeX')


%% XY, plotting
clc;clear all

%data = load('XY/TEstdErhoXY_T_0.750-1.060_32_L-128_PERIODIC.tsv');
data32 = load('XY/TEstdErhoXY_T_0.010-1.900_64_L-32_PERIODIC.tsv');
data8 = load('XY/TEstdErhoXY_T_0.010-1.900_64_L-8_PERIODIC.tsv');


T32 = data32(:,1);
E32 = data32(:,2);
stdE32 = data32(:,3);
rhox32 = data32(:,4); rhoy32 = data32(:,5);

T8 = data8(:,1);
E8 = data8(:,2);
stdE8 = data8(:,3);
rhox8 = data8(:,4); rhoy8 = data8(:,5);



figure(1), clf
subplot(1,2,1)
%plot(T,stdE.^2./T.^2,'.')
plot(T8,E8/8^2, 'b.',T32,E32/32^2, 'r+')
l=legend('$L=8$','$L=32$','location','SouthEast');
set(l,'interpreter','Latex');

set(gca, 'xScale', 'lin', 'fontsize',13);
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$\langle E\rangle/(JN)$', 'interpreter', 'LaTeX');


%
subplot(1,2,2)
CV8=stdE8.^2./T8.^2;
CV32=stdE32.^2./T32.^2;
plot(T8,CV8/max(CV8),'.b',T32,CV32/max(CV32),'+r')
l=legend('$L=8$','$L=32$','location','South');
set(l,'interpreter','Latex');

set(gca, 'yScale', 'lin', 'fontsize',13)%, 'ylim',[0,1])
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$\hat{C}_V$', 'interpreter', 'LaTeX');


%%
figure(2), clf
plot(T8,rhox8, '.b', T8,rhoy8,'ob'), hold on
plot(T32,rhox32, '+r', T32,rhoy32,'xr')
tt=linspace(0,2);
plot(tt, tt*2/pi, 'k')


set(gca, 'yScale', 'lin', 'ylim',[0,1], 'fontsize',13)
grid on;
xlabel('$T/J$', 'interpreter', 'LaTeX');
ylabel('$\langle\rho_{\rm s}\rangle$', 'interpreter', 'LaTeX');

l=legend('$\langle\rho_{\rm s}\rangle_x$, $L=8$',...
         '$\langle\rho_{\rm s}\rangle_y$, $L=8$',...
         '$\langle\rho_{\rm s}\rangle_x$, $L=32$',...
         '$\langle\rho_{\rm s}\rangle_y$, $L=32$',...
         '$2T/(J\pi)$');
set(l,'interpreter','LaTeX');



%% XY finding T_KT
clc, clear all,clf

L=2.^(4:7).';
rhos=zeros(32,length(L));

for i=1:length(L)
    data = load(sprintf('XY/TEstdErhoXY_T_0.750-1.060_32_L-%d_PERIODIC.tsv',L(i)));
    rhos(:,i)=sum(data(:,4:5),2)/2;
end
T = data(:,1);
clear data

t=linspace(0.8,1.2);
c=2/pi;

figure(3)
plot(T,rhos(:,1),'.b',T,rhos(:,2),'*r',...
     T,rhos(:,3),'oc',T,rhos(:,4),'dm')
hold on
plot(t,t*c,'-k')
axis([.75,1.1,.55,.8])

%%
index=16:27;
y=rhos(index, :);
x=T(index);

A=[x, ones(size(x))]\y;
Y=[x, ones(size(x))]*A;

plot(x,Y(:,1),':b',x,Y(:,2),'--r',...
     x,Y(:,3),'-c',x,Y(:,4),'-.m')

l=legend('$L=12$','$L=16$','$L=20$','$L=24$',...
         '$y=2T/\pi$','Fit $L=12$','Fit $L=16$','Fit $L=20$','Fit $L=24$');
set(l,'interpreter','LaTeX')


intercept = ( A(2,:)./(c-A(1,:)) ).';
%%
clc
figure(4),clf,hold on
xx=log(L).^(-2);
plot(xx,intercept, '*')



B1=[xx(1:end), ones(4,1)]\intercept(1:end);
B2=[xx(2:end), ones(3,1)]\intercept(2:end);
B3=[xx(1:3), ones(3,1)]\intercept(1:3);

fprintf('T_KT = %0.3f\n', B3(2));
XX=linspace(0,.17,8).';
AA=[XX,ones(size(XX))];
plot(XX,AA*B1,XX,AA*B2,XX,AA*B3)





%% XY finding T_KT
clc, clear all,clf

L=[20:20:120, 200].';
rhos=zeros(32,length(L));

for i=1:length(L)
    data = load(sprintf('TEstdErhoXY_T_0.950-1.105_32_L-%d_PERIODIC.tsv',L(i)));
    rhos(:,i)=sum(data(:,4:5),2)/2;
end
T = data(:,1);
clear data

t=linspace(0.95,1.1);
c=2/pi;

plot(T,rhos,'.-')
hold on
plot(t,t*c)
%axis([.95,1.1,.6,.7])


%%
clc
index=1:32;
y=rhos(index, :);
x=T(index);



A=[x, ones(size(x))]\y;

intercept = ( A(2,:)./(c-A(1,:)) ).';

clf
xx=log(L).^(-2);
plot(xx,intercept, '*'), hold on

runs=1:6;
B=[xx(runs), ones(size(xx(runs)))]\intercept(runs);

xxx=linspace(0,.12);
plot(xxx, B(1)*xxx+B(2))


fprintf('T_KT = %0.3f\n', B(2));









