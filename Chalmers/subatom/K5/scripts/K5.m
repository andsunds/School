%% intro
clc;clf;clearvars

u=931.504; %MeV

m_137La=136.906350; %u
m_137Ba=136.905815; %u
m_137Cs=136.907074; %u
m_137Xe=136.911742; %u

Q=(m_137Cs-m_137Ba)*u

%% calibration
clc;clf;clearvars

filnamn=cell(3,1);
filnamn{1}='background_97h.txt';
filnamn{2}='207bi_46.5h.txt';
filnamn{3}='137cs_117h.txt';

fil=3;

data=load(filnamn{fil}, '-ascii');
ch=data(:,1);
N=data(:,2);

figure(1)
plot(ch,N); hold on
%set(gca, 'YScale', 'log')

rangemin=[3540, 3710];
rangemax=[3600, 3790];
E=[758;800];%[1]

l=length(E);

x=zeros(l,1);
y=zeros(l,1);

for i=1:l
    [x(i),~, ~]=fitpeak(rangemin(i), rangemax(i), ch,N);
    plot([x(i),x(i)], [1,4e4], 'r')
end

k=[x,ones(size(x))]\E;

figure(2)
e=k(1)*ch+k(2);
plot(e,N)




%[1] M. Röv, 'Fake values', Jour. of Bull Shit, jan. 2222