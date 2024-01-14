RR=importdata("w_right_comparison_to_110.xls");
LR=importdata("w_left_comparison_to_110.xls");
temp1=reshape(RR.data.',1,[])
temp2=reshape(LR.data.',1,[])
A=[temp1,temp2]
S=std(temp1)
S1=std(temp2)
f1=figure()
h= histogram(A)
h.NumBins = 20
xlim([0.60 1])
ylim([0 85])
title(" R values of muscles synergies for full muscle set")

f2=figure()

h2= histogram(temp1)
h2.NumBins = 25
xlim([0.64 1])
title(" R values of muscles synergies for Right muscles")



















































