function [rsqure_c, rsqure_m, rsqure]=funR(data,W,C)
%[Rsqrcond, Rsqrmus, Rsqr]=funvaf(data,w,c) calculates the overall frobineus 
% correlation coefficient between the reconstructed DATA,
%where reconstructed data= W x C.
%It also determines R^2 correlation coeff for the reconstruction of each
%muscle (Rsqrmus) and for each condition (Rsqrcond).


[nmus ncond]=size(data);
[nsyn ndum]=size(C);

%Calculate reconstructed values
ReconAll=W*C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Correlation Coefficients%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate correlation coefficient in the reconstruction of muscles 
%activation patterns for each muscle

for m=1:nmus
    x=data(m,:);
    y=ReconAll(m,:);
    [r,p]=corrcoef(x,y);
    pvalue_m(m)=p(1,2);
    if pvalue_m(m)<0.05
        rsqure_m(m)=r(1,2)^2;
    else
        rsqure_m(m)=NaN;
    end
end

%Calculate correlation coefficient in the reconstruction of muscles 
%activation patterns for each direction

for j=1:ncond
    x=data(:,j);
    y=ReconAll(:,j);
    [r,p]=corrcoef(x,y);
    pvalue_c(j)=p(1,2);
    if pvalue_c(j)<0.05
        rsqure_c(j)=r(1,2)^2;
    else
        rsqure_c(j)=NaN;
    end
end

%Calculate overall R^2 (1x1)

x=data;
y=ReconAll;
[r,p]=corrcoef(x,y);
pvalue=p(1,2);
if pvalue<0.05
    rsqure=r(1,2)^2;
else
    rsqure=NaN;
end



    
