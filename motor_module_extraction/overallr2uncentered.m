function UR=overallr2uncentered(data,ReconData)
%overallr2uncentered
%Calculates overall variability(1x1)

X=cat(3,data,ReconData);    
UR=(sum(sum(prod(X,3))))^2/(sum(sum(data.^2))*sum(sum(ReconData.^2)));
UR=100*UR;