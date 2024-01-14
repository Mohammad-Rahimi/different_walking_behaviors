%this is left and right comparison bar plots:
%% we compared Cs of R and L leg for the same solution Reported Values to R values
%symmetric for right leg
% C(speed,asym,stepwidth)
%speed 0.8  to  1.2  to  1.45
%asym   0   to  15%  to   30#
%width normal  wider     wider
for j=1:2
    if j==1
    data1=importdata("CsMeanMaxRCoA_finalRight.xls");
    side="Right"
    elseif j==2
    data1=importdata("CsMeanMaxRCoA_finalLeft.xls");
    side="Left"    
    end
 list=["C1average","C2average","C2average","C3average","C1max","C2max","C3max",...
     "C4max","C1R","C2R","C3R","C4R","C1CoA","C2CoA","C3CoA","C4CoA","VAF"]
for jj=1:length(list)

%%%%%%
% C(speed,asym,stepwidth)
%% VAF of Cs
for i=1:3
VAFr010(1,i)=(data1.data(i,jj));
VAFr010(2,i)=(data1.data(i+9,jj));
VAFr010(3,i)=(data1.data(i+18,jj));
%%%
VAFr020(1,i)=(data1.data(i+3,jj));
VAFr020(2,i)=(data1.data(i+12,jj));
VAFr020(3,i)=(data1.data(i+21,jj));
%%%
VAFr030(1,i)=(data1.data(i+6,jj));
VAFr030(2,i)=(data1.data(i+15,jj));
VAFr030(3,i)=(data1.data(i+24,jj));
end
avgR(:,1,:)=VAFr010;
avgR(:,2,:)=VAFr020;
avgR(:,3,:)=VAFr030;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coltick={"0.8m/s","1.2m/s","1.45m/s"};
rowtick={"normal","wide","wider"};

% errhigh=maxR-avgR;
% errlow=avgR-minR;

yy = ["0%" "15%" "30%"]

f(2*jj+j)=figure()
%%
ax1 = nexttile;
xx1(:,:)=avgR(1,:,:);


b1=bar3(ax1,xx1,1,'grouped')
view(-20,20);
b1(1).FaceColor='#0D4A70'
b1(2).FaceColor='#226E9C'
b1(3).FaceColor='#3C93C2'

b1(1).LineWidth=1
b1(2).LineWidth=1
b1(3).LineWidth=1

minbar=min(data1.data(:,jj))
maxbar=max(data1.data(:,jj))
zlim([minbar*0.95 min(maxbar*1.05,1)]);
set(gca,'yticklabel',yy)
title("0.8m/s")

%%
ax2 = nexttile;
xx2(:,:)=avgR(2,:,:)


b2=bar3(ax2,xx2,1,'grouped')
view(-20,20);
b2(1).FaceColor='#0D4A70'
b2(2).FaceColor='#226E9C'
b2(3).FaceColor='#3C93C2'
b2(1).LineWidth=1
b2(2).LineWidth=1
b2(3).LineWidth=1


zlim([minbar*0.95 min(maxbar*1.05,1)]);
set(gca,'yticklabel',yy)
title("1.1m/s")

ax3 = nexttile;
xx3(:,:)=avgR(3,:,:)




b3=bar3(ax3,xx3,1,'grouped')
view(-20,20);
b3(1).FaceColor='#0D4A70'
b3(2).FaceColor='#226E9C'
b3(3).FaceColor='#3C93C2'
b3(1).LineWidth=1
b3(2).LineWidth=1
b3(3).LineWidth=1


zlim([minbar*0.95 min(maxbar*1.05,1)]);
set(gca,'yticklabel',yy)
title("1.45m/s")
legend({'normal', 'wide', 'wider'});
sgtitle(strcat(list(jj),side),'Color','red')
set(sgtitle, 'FontSize', 20)

exportgraphics(f(2*jj+j),strcat(list(jj),side,'.jpg'),'Resolution',500)
print(f(2*jj+j), strcat(list(jj),side),'-dpdf','-painters','-bestfit');
savefig(strcat(list(jj),side))
end



end
