%this is left and right comparison bar plots:
%% we compared R of Ws for R and L leg for the same solution Reported Values to R values
%symmetric for right leg
% C(speed,asym,stepwidth)
%speed 0.8  to  1.2  to  1.45
%asym   0   to  15%  to   30#
%width normal  wider     wider
data1=importdata("W_left_to_right_comparison_all_solutions.xls");

%%%%%%
% C(speed,asym,stepwidth)
for i=1:3
VAFr010(1,i)=nanmean(data1.data(i,:));
VAFr010(2,i)=nanmean(data1.data(i+9,:));
VAFr010(3,i)=nanmean(data1.data(i+18,:));
%%%
VAFr020(1,i)=nanmean(data1.data(i+3,:));
VAFr020(2,i)=nanmean(data1.data(i+12,:));
VAFr020(3,i)=nanmean(data1.data(i+21,:));
%%%
VAFr030(1,i)=nanmean(data1.data(i+6,:));
VAFr030(2,i)=nanmean(data1.data(i+15,:));
VAFr030(3,i)=nanmean(data1.data(i+24,:));
end
avgR(:,1,:)=VAFr010;
avgR(:,2,:)=VAFr020;
avgR(:,3,:)=VAFr030;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% min values
for i=1:3
min1(1,i)=nanmin(data1.data(i,:));
min1(2,i)=nanmin(data1.data(i+9,:));
min1(3,i)=nanmin(data1.data(i+18,:));
%%%
min2(1,i)=nanmin(data1.data(i+3,:));
min2(2,i)=nanmin(data1.data(i+12,:));
min2(3,i)=nanmin(data1.data(i+21,:));
%%%
min3(1,i)=nanmin(data1.data(i+6,:));
min3(2,i)=nanmin(data1.data(i+15,:));
min3(3,i)=nanmin(data1.data(i+24,:));
end
minR(:,1,:)=min1;
minR(:,2,:)=min2;
minR(:,3,:)=min3;
%%
%% max values
for i=1:3
max1(1,i)=nanmax(data1.data(i,:));
max1(2,i)=nanmax(data1.data(i+9,:));
max1(3,i)=nanmax(data1.data(i+18,:));
%%%
max2(1,i)=nanmax(data1.data(i+3,:));
max2(2,i)=nanmax(data1.data(i+12,:));
max2(3,i)=nanmax(data1.data(i+21,:));
%%%
max3(1,i)=nanmax(data1.data(i+6,:));
max3(2,i)=nanmax(data1.data(i+15,:));
max3(3,i)=nanmax(data1.data(i+24,:));
end
maxR(:,1,:)=max1;
maxR(:,2,:)=max2;
maxR(:,3,:)=max3;

coltick={"0.8m/s","1.2m/s","1.45m/s"};
rowtick={"normal","wide","wider"};

errhigh=maxR-avgR;
errlow=avgR-minR;

yy = ["0%" "15%" "30%"]

f=figure('Renderer', 'painters', 'Position', [10 10 1600 900])
%%
ax1 = nexttile;
xx1(:,:)=avgR(1,:,:);
xx11(:,:)=minR(1,:,:);
xx111(:,:)=maxR(1,:,:);

b111=bar3(ax1,xx111,1,'grouped')
view(-20,20);
b111(1).FaceColor='#226E9C'
b111(2).FaceColor='#3C93C2'
b111(3).FaceColor='#6CB0D6'
b111(1).LineStyle=':'
b111(2).LineStyle=':'
b111(3).LineStyle=':'
b111(1).LineWidth=1.2
b111(2).LineWidth=1.2
b111(3).LineWidth=1.2
hold on

b1=bar3(ax1,xx1,1,'grouped')
view(-20,20);
b1(1).FaceColor='#226E9C'
b1(2).FaceColor='#3C93C2'
b1(3).FaceColor='#6CB0D6'
b1(1).LineWidth=1
b1(2).LineWidth=1
b1(3).LineWidth=1
hold on


b11=bar3(ax1,xx11,1,'grouped')

minbar=min(min(min(minR)))
zlim([minbar*0.95 1])
set(gca,'yticklabel',yy)
title("0.8m/s")
view(-20,20);
b11(1).FaceColor='#226E9C'
b11(2).FaceColor='#3C93C2'
b11(3).FaceColor='#6CB0D6'
b11(1).LineStyle=':'
b11(2).LineStyle=':'
b11(3).LineStyle=':'
b11(1).LineWidth=1.2
b11(2).LineWidth=1.2
b11(3).LineWidth=1.2
hold off
%%
ax2 = nexttile;
xx2(:,:)=avgR(2,:,:)
xx22(:,:)=minR(2,:,:);
xx222(:,:)=maxR(2,:,:);

b222=bar3(ax2,xx222,1,'grouped')
view(-20,20);
b222(1).FaceColor='#226E9C'
b222(2).FaceColor='#3C93C2'
b222(3).FaceColor='#6CB0D6'
b222(1).LineStyle=':'
b222(2).LineStyle=':'
b222(3).LineStyle=':'
b222(1).LineWidth=1.2
b222(2).LineWidth=1.2
b222(3).LineWidth=1.2
hold on

b2=bar3(ax2,xx2,1,'grouped')
view(-20,20);
b2(1).FaceColor='#226E9C'
b2(2).FaceColor='#3C93C2'
b2(3).FaceColor='#6CB0D6'
b2(1).LineWidth=1
b2(2).LineWidth=1
b2(3).LineWidth=1
hold on

b22=bar3(ax2,xx22,1,'grouped')
view(-20,20);
b22(1).FaceColor='#226E9C'
b22(2).FaceColor='#3C93C2'
b22(3).FaceColor='#6CB0D6'
b22(1).LineStyle=':'
b22(2).LineStyle=':'
b22(3).LineStyle=':'
b22(1).LineWidth=1.2
b22(2).LineWidth=1.2
b22(3).LineWidth=1.2

zlim([minbar*0.95 1])
set(gca,'yticklabel',yy)
title("1.2m/s")

hold off
ax3 = nexttile;
xx3(:,:)=avgR(3,:,:)
xx33(:,:)=minR(3,:,:);
xx333(:,:)=maxR(3,:,:);

b333=bar3(ax3,xx333,1,'grouped')
view(-20,20);
b333(1).FaceColor='#226E9C'
b333(2).FaceColor='#3C93C2'
b333(3).FaceColor='#6CB0D6'
b333(1).LineStyle=':'
b333(2).LineStyle=':'
b333(3).LineStyle=':'
b333(1).LineWidth=1.2
b333(2).LineWidth=1.2
b333(3).LineWidth=1.2
hold on


b3=bar3(ax3,xx3,1,'grouped')
view(-20,20);
b3(1).FaceColor='#226E9C'
b3(2).FaceColor='#3C93C2'
b3(3).FaceColor='#6CB0D6'
b3(1).LineWidth=1
b3(2).LineWidth=1
b3(3).LineWidth=1
hold on
b33=bar3(ax3,xx33,1,'grouped')
view(-20,20);
b33(1).FaceColor='#226E9C'
b33(2).FaceColor='#3C93C2'
b33(3).FaceColor='#6CB0D6'
b33(1).LineStyle=':'
b33(2).LineStyle=':'
b33(3).LineStyle=':'
b33(1).LineWidth=1.2
b33(2).LineWidth=1.2
b33(3).LineWidth=1.2

zlim([minbar*0.95 1])
set(gca,'yticklabel',yy)
title("1.45m/s")
legend({'normal', 'wide', 'wider'});
sgtitle('R of right to right of 110 ','Color','red')
set(sgtitle, 'FontSize', 20)






