%% Kalibrering med Na-topparna
clc;clf;clear all

addpath('Data/')

data=load('kalibrering4.tsv');

lambda=data(:,1);
I=data(:,2);

plot(lambda,I)


L=length(I);

[M, i]=max(I(1:floor(L/2)));

fprintf('First Na-peak at lambda = %3.5f. \n\n', lambda(i))

%% Kalibrering med Na-topparna
clc;clf;clear all
addpath('Data/obsolet')

data1=load('test_002.lvm');
data2=load('test_003.lvm');

lambda1=data1(:,2);
I1=data1(:,1);
lambda2=data2(:,2);
I2=data2(:,1);

plot(lambda1,I1,lambda2,I2 )

lambda=mean([data1(:,2), data2(:,2)]');
I=mean([data1(:,1), data2(:,1)]');

figure(2), clf
plot(lambda,I)


%% Bakgrundsmätning, halogen
clc;clf;clear all


a=dir('Data/halogen_wasser_kompl*');

L=length(a);

X=zeros(211,L);
Y=zeros(211,L);

for i=1:L
    file=a(i).name;
    data=load(file);
    X(:,i)=data(:,2);
    Y(:,i)=data(:,1);
end
la=mean(X,2);
Ia=mean(Y,2);

%plot(la,Ia)

b=dir('Data/halogen_wasser_0*');

L=length(b);

X=zeros(701,L);
Y=zeros(701,L);

for i=1:L
    file=b(i).name;
    data=load(file);
    X(:,i)=data(:,2);
    Y(:,i)=data(:,1);
end

lb=mean(X,2);
Ib=mean(Y,2);


L_background=[la(1:200); (la(201)+lb(1:10))/2; lb(11:end)];
I_background=[Ia(1:200); (Ia(201)+Ib(1:10))/2; Ib(11:end)];

%Av eller på för absorption eller emission
I_background = (I_background-min(I_background));

figure(1);clf
plot(L_background,I_background, '-k')


load('constants.mat', 'h', 'c', 'k_B')


B_lambda = @(lambda, T) 2*h*c^2.*lambda.^(-5) ./( exp(h*c./(lambda*k_B.*T)) - 1 );

hold on

X=2600:200:3200;
l=linspace(200,1200);
for T=X;
y=B_lambda(l*1e-9, T);

%plot(l, y/max(y)*max(I_background), '-.')
plot(l, y/2.25e12, '-.')
end
hold off



%Den här är bra att använda om man vill plotta flera grafter
%där ett värde har ändrats för varje graf.
legendCell = regexp(sprintf('T=%d',X),...
                    'T=-?+\d+\.?\d*', 'match');
                    
legendCell = regexprep(legendCell,'\.','\{,\}'); %ändrar till decimalKOMMA
legendCell = strcat('$',legendCell, '\,\mathrm{K}$'); %'$' i början och slutet för att få LaTeX-kod

l1=legend([{'M\"a{}tdata'},legendCell]);
set(l1, 'interpreter', 'latex', 'fontsize', 15, 'location', 'best')
set(gca,'fontsize', 15, 'yscale', 'log', 'xscale', 'lin')

ylabel('Intensitet /[godt. enhet]', 'interpreter', 'latex', 'fontsize', 20)
xlabel('$\lambda$/[nm]', 'interpreter', 'latex', 'fontsize', 20)

axis([200, 1200, 1e-3, 1])


figure(2);
clf;clc

L_min=340;

L=L_background(L_background>L_min);
I=I_background(L_background>L_min);


T=3200;
Y=B_lambda(L*1e-9, T);
Y=Y/max(Y);

Z=-log(abs(I./Y));


index=find( (Z>-7) );

plot(L,exp(Z), '--', 'linewidth', 3)
axis([L_min, 1100, 1, 3e2])

hold on
S=linspace(L_min, 1100, 10000);
% deg_abs=34, deg_em=17
P_korr=polyfit(L(index), Z(index), 17); %Verklig I = I_mätt.*exp(P_korr)

h=plot(S, exp(polyval(P_korr, S)) );
%title('Korrektion')

% 
% for i=30:50
%     i
%     P_korr=polyfit(L(index), Z(index), i);
%     set(h, 'Ydata', exp(polyval(P_korr, S)))
%     pause()
% end
%
set(gca, 'fontsize', 15, 'yscale', 'log')
l2=legend('Uppm\"a{}tt data', 'Anpassad korrektionsfaktor', 'location', 'best');
set(l2, 'interpreter', 'LaTeX', 'fontsize', 13)

ylabel('$I_\mathrm{Planck}/I_\mathrm{spektrometer}$', 'interpreter', 'latex', 'fontsize', 20)
xlabel('$\lambda$/[nm]', 'interpreter', 'latex', 'fontsize', 20)

figure(3);clf
hold on
plot(L,(abs(I)), '--')
plot(L,(Y), '-.', 'linewidth', 3)

plot(L,...
    (abs(I)) .* exp(polyval(P_korr, L)), '-k','linewidth', 2)

axis([L_min, 1100, 1e-2, 1.2e0])
%title('Korrigerad int.')
set(gca, 'fontsize', 14, 'yscale', 'log')

l3=legend('Uppm\"a{}tt data', 'Plancks str\aa{}lningslag','Korrigerad data', 'location', 'best');
set(l3, 'interpreter', 'LaTeX', 'fontsize', 13)

ylabel('Intensitet /[godt. enhet]', 'interpreter', 'latex', 'fontsize', 20)
xlabel('$\lambda$/[nm]', 'interpreter', 'latex', 'fontsize', 20)


%%
clc
%L_background=L;
%I_background=I;

%save('Data/halogen_backgound_absorption.mat','L_background','I_background', 'L_min', 'P_korr')
%save('Data/halogen_backgound_emission.mat','L_background','I_background', 'L_min', 'P_korr')

disp('Save successfull!')












