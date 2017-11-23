%%
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

V_cell=cell(N);
I_cell=cell(N);

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

%% DEBUG
GPV=spec{i}.gaussParVel
npk=length(GPV)/3;
V=GPV( 3*(1:npk)-1 );
Int=GPV( 3*(1:npk)-2 );
plot(V,Int,'ko')

%% DEBUG
x=linspace(-200,200);
y1=GPV(1)*exp(-((x-GPV(2))/GPV(3)).^2/2);
y2=GPV(4)*exp(-((x-GPV(5))/GPV(6)).^2/2);

plot(x,y1, x,y2,x,y1+y2)

%%
clc;clf

V0=220; %km/s (project guide on SALSA web page)
R0=8.5; %kpc  (project guide on SALSA web page)

V=zeros(N,1);
R_min=zeros(N,1);

for i = 1:N
    R_min(i)=R0*sind(longspec(i));
    V_max=max(V_cell{i});
    V(i)=V_max+V0*sind(longspec(i));
    if longspec(i)>90
        V(i)=NaN;
    end
end


plot(R_min,V,'o')
axis([0,8.5,0,250])




















