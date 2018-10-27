function [ df ] = hyp_odefun(f, tau, alpha, dx)

ff=reshape(f, [],2);

df=zeros(size(ff));

df(:,1)=ff(:,2);

df(:,2)=(alpha*[0; diff(ff(:,1),2); 0]/dx-ff(:,1))/tau;

df=reshape(df, [],1);


end

