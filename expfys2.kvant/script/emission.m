%% Hg-plot, egna mätningar
clc;clear all;figure(1), clf

addpath('Data/')
%load('Data/halogen_backgound_20151213.mat', 'P_korr', 'L_min')
load('Data/halogen_backgound_emission.mat', 'P_korr', 'L_min')



threshold1=6e-4;
threshold2=.9e-3;

hold on
plot([L_min, 1100], threshold1*[1,1], '-.k',...
     [L_min, 1100], threshold2*[1,1], '--k')

%Visar korrektionen
x=linspace(L_min, 1100,1000);
y=threshold2*exp(polyval(P_korr, x));
plot(x,y, ':k')

%P_korr=0;
a=dir('Data/Hg_20151212_00*.lvm');
bra=[2];
legendCell1=cell(size(bra));
normfaktor=0;
for i=1:length(bra)
data1=load(a(bra(i)).name);

L1=data1(:,1);
I1=data1(:,2);

index1=find(I1>threshold2);
I1(index1)=I1(index1).*exp(polyval(P_korr, L1(index1)));
plot(L1,I1)

legendCell1{i}=num2str(bra(i));
normfaktor=normfaktor+max(I1(L1<500))/length(bra);
end

axis([L_min, 1050, 4e-4, 10])
grid on
set(gca, 'fontsize', 15, 'yscale', 'log')

%

b=dir('Data/NIST_persistent_lines_Hg.tsv');

legendCell2=cell(size(b'));

for i=1:length(b)
    data=load(b(i).name);
    I2=data(:,1)/1000; I2=I2/max(I2);
    L2=data(:,2)/10;% Å -> nm
    plot(L2, I2*normfaktor, 'o', 'markersize', 10)
    
    r=regexp(b(i).name,'[\_\.]' ,'split');
    legendCell2{i}=[r{end-1}, ' "persistent lines"'];
    %legendCell2{i}='"persistent lines"';
end


data=load('Data/NIST_strong_lines_Hg.tsv');
I=data(:,1); I=I/max(I);
L3=data(:,2)/10;% Å -> nm
plot(L3, I*normfaktor, '*k')


toppar=find( ([0; diff(I1)]>0) & ([diff(I1); 0]<0) & (I1>threshold2*1.5));
plot(L1(toppar), I1(toppar), 'vb')



l=legend([{'Tr\"o{}skel 1','Tr\"o{}skel 2', 'Korrektion'},...
       {'Spektrum'},...%legendCell1, ...
       legendCell2,...
       {'Hg "strong lines"', 'Detekterade toppar'}],...
       'location', 'northeast');

set(l, 'interpreter', 'latex', 'fontsize', 13)


ylabel('Intensitet /[godt. enhet]', 'interpreter', 'latex', 'fontsize', 20)
xlabel('$\lambda$/[nm]', 'interpreter', 'latex', 'fontsize', 20)

%% Energier
clc;
%figure(2),
clf; 
hold on
load('constants.mat', 'h', 'c', 'e')

E1=h*c./(L1(toppar)*1e-9)/e;
E3=h*c./(L3*1e-9)/e;

A=[L1(toppar), E1,   1e-5./(L1(toppar)*1e-9)];

%repmat(E1, 1, length(E1)+1) - [zeros(size(E1)), E1.'.*ones(length(E1))]
% v- nm,     v- eV,     v- 10^3 cm^-1
str1=sprintf('%3.2f  \t& %3.3f \t& %3.2f \\\\ \n',A.')


save('toppar.mat', 'A', '-mat')


xline=[-1 1];
for i=1:length(E1)
plot(xline, E1(i)*ones(size(xline)), '.-k')
end

for i=1:length(E3)
plot(xline+ 0*(i+2-length(E1)), E3(i)*ones(size(xline)), '.r')
end

% yyaxis left
set(gca, 'xtick', [], 'ytick', 1:.1:3.5, 'ylim', [0,3.5])

% yyaxis right
% set(gca, 'ydir', 'reverse','ytick', e*1e9/h/c./fliplr([0:0.5:3.5]),'ylim', [200, 1100])
%%
clc;clf;clear all

load('constants.mat', 'h', 'c', 'e');

data=load('Data/NIST_energylevels_Hg.tsv');
[data-data(1),data]

%% Kr-plot, egna mätningar
clc;clear all; clf

addpath('Data/obsolet/')
%load('Data/halogen_backgound_20151213.mat', 'P_korr', 'L_min')
load('Data/halogen_backgound_emission.mat', 'P_korr', 'L_min')



threshold1=6e-4;
threshold2=7e-3;

hold on
plot([L_min, 1100], threshold2*[1,1], '-.k',...
     [L_min, 1100], threshold1*[1,1], '--k')

%Visar korrektionen
x=linspace(L_min, 1100,1000);
y=threshold2*exp(polyval(P_korr, x));
plot(x,y, ':k')

%P_korr=0;
a=dir('Data/obsolet/data_Kr_0.tsv');
bra=[1];
legendCell1=cell(size(bra));
normfaktor=0;
for i=1:length(bra)
data1=load(a(bra(i)).name);

L1=data1(:,1);
I1=data1(:,2);

index1=find(I1>threshold2);
I1(index1)=I1(index1).*exp(polyval(P_korr, L1(index1)));
plot(L1,I1)

legendCell1{i}=num2str(bra(i));
normfaktor=normfaktor+max(I1(L1>500))/length(bra);
end

axis([L_min, 1050, 4e-4, 10])
grid on
set(gca, 'fontsize', 15, 'yscale', 'log')

%

b=dir('Data/NIST_persistent_lines_Kr.tsv');

legendCell2=cell(size(b'));

for i=1:length(b)
    data=load(b(i).name);
    I2=data(:,1)/1000; I2=I2/max(I2);
    L2=data(:,2)/10;% Å -> nm
    plot(L2, I2*normfaktor, 'o', 'markersize', 10)
    
    r=regexp(b(i).name,'[\_\.]' ,'split');
    legendCell2{i}=[r{end-1}, ' "persistent"'];
    %legendCell2{i}='"persistent lines"';
end


data=load('Data/NIST_strong_lines_Kr.tsv');
I=data(:,1); I=I/max(I);
L3=data(:,2)/10;% Å -> nm
plot(L3, I*normfaktor, '*k')

l=legend([{'tr\"o{}skel 1','tr\"o{}skel 2', 'Korrektion'},...
       legendCell1, legendCell2,...
       ],...
       'location', 'northeast');

set(l, 'interpreter', 'latex', 'fontsize', 14)




