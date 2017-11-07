%%
clc;clf;clear
direc = './2017Nov07/';
files = dir ([direc, '*.fits']);

for i = 1:length( files )
    fname = [ direc, files(i).name ];
    spec{i}= SalsaSpectrum(fname);
    longspec(i) = spec{i}.getKeyword ('CRVAL2');
    l = longspec(i);
    if(l >0 && l <40)
        spec{i}. fitBaseline ([ -240,-205 -140 -100 120 220] , 'vel' ,3);
    elseif(l >=40 && l <90)
        spec{i}.fitBaseline([-230,-200 -150 -110 100 220] , 'vel' ,3);
    elseif (l >=90 && l <180)
        spec{i}.fitBaseline([ -230,-190 -130 -110 50 220] , 'vel' ,3);
    elseif (l>=180)
        spec{i}.fitBaseline([ -230,-180 -120 -40 60 220] , 'vel' ,3);
    end

spec{i}.showBaseline()
spec{i}.subtractBaseline()
spec{i}.fitGaussians()
spec{i}.plot()
GPV=spec{i}.gaussParVel;
npk=length(GPV)/3;
V=GPV( 3*(1:npk)-1 );
Int=GPV( 3*(1:npk)-2 );
plot(V,Int,'ko')
pause()
end % end loop over spectra

%%
GPV=spec{i}.gaussParVel
npk=length(GPV)/3;
V=GPV( 3*(1:npk)-1 );
Int=GPV( 3*(1:npk)-2 );
plot(V,Int,'ko')

%%
x=linspace(-200,200);
y1=GPV(1)*exp(-((x-GPV(2))/GPV(3)).^2/2);
y2=GPV(4)*exp(-((x-GPV(5))/GPV(6)).^2/2);

plot(x,y1, x,y2,x,y1+y2)
