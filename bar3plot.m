%number of synergies for left leg across different speeds and stepwidth
%symmetric for right leg
% C(speed,asym,stepwidth)
%speed 0.8  to  1.1  to  1.45
%asym   0   to  15%  to   30#
%width normal  wider     wider
%% i manually made these two files from this file :compareCs_costs2 in the format of synnumbers file
data=importdata("recons_VAF_l.xls");
data1=importdata("recons_VAF_r.xlsx")
% zr010=[4,4,4;4,4,5;4,4,5];
% zr020=[4,4,5;4,4,5;4,4,6];
% zr030=[4,4,5;4,4,6;4,4,6];
% numsynr(:,1,:)=zr010;
% numsynr(:,2,:)=zr020;
% numsynr(:,3,:)=zr030;
% 
% VAFr010=[98.0686,97.8497,97.7689;97.8692,97.7173,97.0411;97.6617,97.5493,96.8212];
% VAFr020=[97.9368,97.8243,96.8398;97.8845,97.7866,96.3378;97.7139,97.701,95.8327];
% VAFr030=[98.0021,97.8271,96.6194;98.0265,97.7819,96.0249;97.6728,97.5991,95.5077];
% VAFr(:,1,:)=VAFr010;
% VAFr(:,2,:)=VAFr020;
% VAFr(:,3,:)=VAFr030;
%%%%%%

zl010(1,:)=data.data(1,1:3);
zl010(2,:)=data.data(1,10:12);
zl010(3,:)=data.data(1,19:21);
%%%%
zl020(1,:)=data.data(1,4:6);
zl020(2,:)=data.data(1,13:15);
zl020(3,:)=data.data(1,22:24);
%%%%;
zl030(1,:)=data.data(1,7:9);
zl030(2,:)=data.data(1,16:18);
zl030(3,:)=data.data(1,25:27);

numsynl(:,1,:)=zl010;
numsynl(:,2,:)=zl020;
numsynl(:,3,:)=zl030;


VAFl010(1,:)=data.data(1,1:3);
VAFl010(2,:)=data.data(1,10:12);
VAFl010(3,:)=data.data(1,19:21);

VAFl020(1,:)=data.data(1,4:6);
VAFl020(2,:)=data.data(1,13:15);
VAFl020(3,:)=data.data(1,22:24);
%%%
VAFl030(1,:)=data.data(1,7:9);
VAFl030(2,:)=data.data(1,16:18);
VAFl030(3,:)=data.data(1,25:27);

VAFl(:,1,:)=VAFl010*100;
VAFl(:,2,:)=VAFl020*100;
VAFl(:,3,:)=VAFl030*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zr010(1,:)=data1.data(1,1:3);
zr010(2,:)=data1.data(1,10:12);
zr010(3,:)=data1.data(1,19:21);
%%%%
zr020(1,:)=data1.data(1,4:6);
zr020(2,:)=data1.data(1,13:15);
zr020(3,:)=data1.data(1,22:24);
%%%%;
zr030(1,:)=data1.data(1,7:9);
zr030(2,:)=data1.data(1,16:18);
zr030(3,:)=data1.data(1,25:27);

numsynr(:,1,:)=zr010;
numsynr(:,2,:)=zr020;
numsynr(:,3,:)=zr030;


VAFr010(1,:)=data1.data(1,1:3);
VAFr010(2,:)=data1.data(1,10:12);
VAFr010(3,:)=data1.data(1,19:21);

VAFr020(1,:)=data1.data(1,4:6);
VAFr020(2,:)=data1.data(1,13:15);
VAFr020(3,:)=data1.data(1,22:24);
%%%
VAFr030(1,:)=data1.data(1,7:9);
VAFr030(2,:)=data1.data(1,16:18);
VAFr030(3,:)=data1.data(1,25:27);

VAFr(:,1,:)=VAFr010*100;
VAFr(:,2,:)=VAFr020*100;
VAFr(:,3,:)=VAFr030*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



coltick={"0.8m/s","1.1m/s","1.45m/s"};
rowtick={"normal","wide","wider"};

%%%%
% [X,Y] = meshgrid(1:3)
% contour3([X,Y],numsynr(X,Y,1))
%  x=[1;2;3]
%  y=[1;2;3]
%  z=[1;2;3]
% [X,Y,Z]=meshgrid(x,y,z)
%  %contour3(x,y,numsynr(:,:,3)); 
%  %surf(x,y,numsynr(:,:,3)); 
%  %pdeplot3D([X,Y,Z],numsynr)
%  contourf(x,y,(numsynr(:,:,3))')
%  colorbar
 %F = scatteredInterpolant(x,y,z,numsynr(:,:,3)) 
yy = ["0%" "15%" "30%"]

f=figure()

ax1 = nexttile;
xx1(:,:)=numsynl(1,:,:)


b1=bar3(ax1,xx1,1,'grouped')
set(gca,'yticklabel',yy)
zlim([0 4])
title("0.8m/s")
b1(1).FaceColor='#8B0000'
b1(2).FaceColor='#FF4500'
b1(3).FaceColor='#FF8C00'

xx2(:,:)=numsynl(2,:,:)
ax2 = nexttile;
b2=bar3(ax2,xx2,1,'grouped')
zlim([0 4])
set(gca,'yticklabel',yy)
title("1.1m/s")
b2(1).FaceColor='#8B0000'
b2(2).FaceColor='#FF4500'
b2(3).FaceColor='#FF8C00'

xx3(:,:)=numsynl(3,:,:)
ax3 = nexttile;
b3=bar3(ax3,xx3,1,'grouped')
zlim([0 4])
set(gca,'yticklabel',yy)
title("1.45m/s")
b3(1).FaceColor='#8B0000'
b3(2).FaceColor='#FF4500'
b3(3).FaceColor='#FF8C00'

legend({'normal', 'wide', 'wider'});
sgtitle('Left Leg','Color','red')
set(sgtitle, 'FontSize', 20)
%%%%
f1=figure()
zlim([0 4])
ax1 = nexttile;
xx1(:,:)=numsynr(1,:,:)
b1=bar3(ax1,xx1,1,'grouped')
zlim([0 4])
set(gca,'yticklabel',yy)
title("0.8m/s")
view(-20,20);
b1(1).FaceColor='#8B0000'
b1(2).FaceColor='#FF4500'
b1(3).FaceColor='#FF8C00'

xx2(:,:)=numsynr(2,:,:)
ax2 = nexttile;
b2=bar3(ax2,xx2,1,'grouped')
zlim([0 4])
set(gca,'yticklabel',yy)
title("1.1m/s")
view(-20,20);
b2(1).FaceColor='#8B0000'
b2(2).FaceColor='#FF4500'
b2(3).FaceColor='#FF8C00'


xx3(:,:)=numsynr(3,:,:)
ax3 = nexttile;
b3= bar3(ax3,xx3,1,'grouped')
zlim([0 4])
set(gca,'yticklabel',yy)
title("1.45m/s")
view(-20,20);
b3(1).FaceColor='#8B0000'
b3(2).FaceColor='#FF4500'
b3(3).FaceColor='#FF8C00'

legend({'normal', 'wide', 'wider'});
sgtitle('Right Leg','Color','red')
set(sgtitle, 'FontSize', 20)
%LimitsX = xlim; LimitsY = ylim;
%sgtitle('Left', 'HorizontalAlignment', 'left', 'position', [10, 10]);
% for i=1:3
%      for j=1:3
%          for k=1:3
%             
% point([i,j,k],numsynr(i,j,k)) 
% hold on 
% box on
%      
%      end
%      end
% end
%contour(X,Y,Z,50)
f2=figure()

ax1 = nexttile;
xx1(:,:)=VAFl(1,:,:)
b1=bar3(ax1,xx1,1,'grouped')
set(gca,'yticklabel',yy)
zlim([88 96])
title("0.8m/s")
b1(1).FaceColor='#8B0000'
b1(2).FaceColor='#FF4500'
b1(3).FaceColor='#FF8C00'

xx2(:,:)=VAFl(2,:,:)
ax2 = nexttile;
b2=bar3(ax2,xx2,1,'grouped')
set(gca,'yticklabel',yy)
zlim([88 96])
title("1.1m/s")
b2(1).FaceColor='#8B0000'
b2(2).FaceColor='#FF4500'
b2(3).FaceColor='#FF8C00'

xx3(:,:)=VAFl(3,:,:)
ax3 = nexttile;
b3=bar3(ax3,xx3,1,'grouped')
set(gca,'yticklabel',yy)
zlim([88 96])
title("1.45m/s")
b3(1).FaceColor='#8B0000'
b3(2).FaceColor='#FF4500'
b3(3).FaceColor='#FF8C00'
legend({'normal', 'wide', 'wider'});
sgtitle('Left Leg','Color','red')
set(sgtitle, 'FontSize', 20)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3=figure()

ax1 = nexttile;
xx1(:,:)=VAFr(1,:,:)
b1=bar3(ax1,xx1,1,'grouped')
set(gca,'yticklabel',yy)
zlim([88 96])
title("0.8m/s")
b1(1).FaceColor='#8B0000'
b1(2).FaceColor='#FF4500'
b1(3).FaceColor='#FF8C00'

xx2(:,:)=VAFr(2,:,:)
ax2 = nexttile;
b2=bar3(ax2,xx2,1,'grouped')
set(gca,'yticklabel',yy)
zlim([88 96])
title("1.1m/s")
b2(1).FaceColor='#8B0000'
b2(2).FaceColor='#FF4500'
b2(3).FaceColor='#FF8C00'

xx3(:,:)=VAFr(3,:,:)
ax3 = nexttile;
b3=bar3(ax3,xx3,1,'grouped')
set(gca,'yticklabel',yy)
zlim([88 96])
title("1.45m/s")
b3(1).FaceColor='#8B0000'
b3(2).FaceColor='#FF4500'
b3(3).FaceColor='#FF8C00'

legend({'normal', 'wide', 'wider'});
sgtitle('Right Leg','Color','red')
set(sgtitle, 'FontSize', 20)

print(f, 'fsub','-dpdf','-painters','-bestfit');
print(f1, 'f1sub','-dpdf','-painters','-bestfit');
print(f2, 'f2sub','-dpdf','-painters','-bestfit');
print(f3, 'f3sub','-dpdf','-painters','-bestfit');

