%% /2017Nov07/ DON'T TOUCH!!!
clc;clf;clear
direc = './2017Nov07/';
files = dir ([direc, '*.fits']);
N=length(files);

%By visual inspection, the automatic fit fails for these indeses
i_bad=[4,6,8,9,12];
%We therefore have to specify the fit paramters manually
%fit params are specified by [height, location, width, ...]
fitparam_manual={[60,-35,1, 47,-6,5],[35,-50,20, 40,-6,10, 10,-100,10]...
    [50,0,10, 30,-60,10, 10,-100,10], [50,-5,10, 25,-50,5, 25,-75,5], [45,8,5, 35,26,5, 20,-60,15]};

V_cell=cell(N,1);
I_cell=cell(N,1);

for i = 1:N
    fname = [ direc, files(i).name ];
    spec{i}= SalsaSpectrum(fname);
    longspec(i) = spec{i}.getKeyword ('CRVAL2');
    l = longspec(i);
    
    if(l >0 && l <40)
        spec{i}.fitBaseline([-240,-205 -140 -100 120 220] , 'vel' ,3);
    elseif(l >=40 && l <90)
        spec{i}.fitBaseline([-230,-200 -150 -110 100 220] , 'vel' ,3);
    elseif (l >=90 && l <180)
        spec{i}.fitBaseline([-230,-190 -130 -110 50 220] , 'vel' ,3);
    elseif (l>=180)
        spec{i}.fitBaseline([-230,-180 -120 -40 60 220] , 'vel' ,3);
    end
    

    %spec{i}.showBaseline()
    spec{i}.subtractBaseline();
    if sum(i_bad==i)
        spec{i}.fitGaussians(fitparam_manual{i_bad==i});
    else
        spec{i}.fitGaussians();
    end
    %spec{i}.plot()
    GPV=spec{i}.gaussParVel;
    npk=length(GPV)/3;
    V=GPV( 3*(1:npk)-1 );  V_cell{i}=V;
    Int=GPV( 3*(1:npk)-2 );I_cell{i}=Int;
    %plot(V,Int,'ko')
    %pause()
end % end loop over spectra

%save([direc,'extracted_data.mat'],'V_cell','I_cell', 'longspec')

%% /2017Nov23/ DON'T TOUCH!!!
clc;clf;clear
direc = './2017Nov23/';
files = dir ([direc, '*.fits']);
N=length(files);

%By visual inspection, the automatic fit fails for these indeses
i_bad=[1,2,9,12];
i_very_bad=[6];
%We therefore have to specify the fit paramters manually
%fit params are specified by [height, location, width, ...]
fitparam_manual={[55,0,10, 25,-45,5, 20,-70,10],[45,0,5, 15,-80,5, 15,-23,10]...
    [50,10,20, 20,-60,20], [45,10,10, 15,-40,5, 20,90,10, 25,80,5]};

V_cell=cell(N,1);
I_cell=cell(N,1);

for i = 1:N
    fname = [ direc, files(i).name ];
    spec{i}= SalsaSpectrum(fname);
    longspec(i) = spec{i}.getKeyword ('CRVAL2');
    l = longspec(i);
    
    if(l >0 && l <40)
        spec{i}.fitBaseline([-240,-205 -140 -100 120 220] , 'vel' ,3);
    elseif(l >=40 && l <90)
        spec{i}.fitBaseline([-230,-200 -150 -110 100 220] , 'vel' ,3);
    elseif (l >=90 && l <180)
        spec{i}.fitBaseline([-230,-190 -130 -110 50 220] , 'vel' ,3);
    elseif (l>=180)
        spec{i}.fitBaseline([-230,-180 -120 -40 60 220] , 'vel' ,3);
    end
    

    %spec{i}.showBaseline()
    spec{i}.subtractBaseline();
    if sum(i_bad==i)
        spec{i}.fitGaussians(fitparam_manual{i_bad==i});
    else
        spec{i}.fitGaussians();
    end
    
    GPV=spec{i}.gaussParVel;
    npk=length(GPV)/3;
    V=GPV( 3*(1:npk)-1 );  V_cell{i}=V;
    Int=GPV( 3*(1:npk)-2 );I_cell{i}=Int;
    
    if sum(i_very_bad==i)
        V=NaN;
        Int=NaN;
        longspec(i)=NaN;
    end
    spec{i}.plot()
    plot(V,Int,'ko')
    %pause()
end % end loop over spectra

%save([direc,'extracted_data.mat'],'V_cell','I_cell', 'longspec')



%% Extracts all data of interest
clc;clf;clear

V_CELL=cell(0);
GLON=[];
N=0;


direc= {'./2017Nov07/','./2017Nov23/'};
for J=1:length(direc)
    load([direc{J},'extracted_data.mat'])
    N=N+length(longspec);
    V_CELL=[V_CELL; V_cell];
    GLON=[GLON,longspec];
end

%save('ALL_extracted_data.mat', 'V_CELL','GLON','N')

%% ROTATION CURVE
clc;clf;clear
load('ALL_extracted_data.mat');
V0=220; %km/s (project guide on SALSA web page)
R0=8.5; %kpc  (project guide on SALSA web page)

V=zeros(N,1);
R_min=zeros(N,1);

for i = 1:N
    R_min(i)=R0*sind(GLON(i));
    V_max=max(V_CELL{i});
    V(i)=V_max+V0*sind(GLON(i));
    if GLON(i)>90
        V(i)=NaN;
    end
end

plot(R_min,V,'o'), hold on
%plot([0,9],[210,210],':k')
axis([0,9,0,250])

xlabel('$R$ /[kpc]', 'interpreter','latex')
ylabel('$V$ /[km/s]', 'interpreter','latex')

set(gca,'fontsize',14)
grid on


%% 
clc;clf;clear

load('ALL_extracted_data.mat');
V0=220; %km/s (project guide on SALSA web page)
R0=8.5; %kpc  (project guide on SALSA web page)

%Initializing as empy vector to then be appended in the for loop, this is
%because we don't know how many cloudes there are.
R=[]; %dist from GC
d=[]; %dist from us
i_bad=[]; %ambigous positions
l=[]; %galactic longitude

for i = 1:N
    Vr=V_CELL{i};
    Ri=R0*V0*sind(GLON(i))./(V0*sind(GLON(i))+Vr); %finding all the R's for this long
    R=[R,Ri]; %appending to the final vecor
    
    %distance from us is di1+-di2
    di1=R0*cosd(GLON(i));
    di2=sqrt(Ri.^2-R0.^2*sind(GLON(i)).^2);
    %d is ambigous if di2<di1
    i_bad_i=find(di2<di1);
    i_bad=[i_bad, length(d)+i_bad_i];
    d=[d, di1+di2];
    %adds and replicates the long's for all clods in this longitude
    l=[l,repmat(GLON(i),1,length(Ri))];
end

% some poits are very bad, we don't want them
i_very_bad=[find(imag(d)>0), find(d<0)];
d(i_very_bad)=NaN;

%Fixed ambiguity on l=50 deg
i_fixed=33;
l_fixed=l(i_fixed);
R_fixed=R(i_fixed);
d_fixed=R0*cosd(l_fixed)-sqrt(R_fixed.^2-R0.^2*sind(l_fixed).^2);
X_fixed=d_fixed.*sind(l_fixed);
Y_fixed=R0-d_fixed.*cosd(l_fixed);


%positions in the galaxy
X=d.*sind(l);
Y=R0-d.*cosd(l);

%plotting clod poitions
plot(X,Y,'o')
hold on
%plotting the ambigous positions
plot(X(i_bad),Y(i_bad),'rx')
plot([X_fixed,X(i_fixed)],[Y_fixed,Y(i_fixed)],'o')

%Drawing our position
plot(0,R0,'k.','markersize',20)

%Drawing circles
t=linspace(-pi,pi);
r1=8.5;
x1=r1*cos(t);y1=r1*sin(t);
r2=12;
x2=r2*cos(t);y2=r2*sin(t);
r3=10;
x3=r3*cos(t);y3=r3*sin(t);
r4=8;
x4=r4*cos(t);y4=r4*sin(t);
r5=14;
x5=r5*cos(t);y5=r5*sin(t);
r6=16;
x6=r6*cos(t);y6=r6*sin(t);
plot(x1,y1,'-.k',x2,y2,'k:',x3,y3,'k:',x4,y4,'k:',x5,y5,'k:',x6,y6,'k:')

%Drawing line of sight for l=50
r=[0,20];
L=l(i_fixed);
x=r*sind(L);
y=R0-r*cosd(L);
plot(x,y,'--k')

axis equal
axis([-5,15, -15,20])
%grid on

xlabel('$X$ /[kpc]', 'interpreter','latex')
ylabel('$Y$ /[kpc]', 'interpreter','latex')

set(gca,'fontsize',14)










