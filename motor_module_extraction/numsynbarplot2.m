%number of synergies for left leg across different speeds and stepwidth
%symmetric for right leg
% C(speed,asym,stepwidth)
%speed 0.8  to  1.1  to  1.45
%asym   0   to  15%  to   30#
%width normal  wider     wider
data=importdata("synnumbers_l.xls");
data1=importdata("synnumbers_r.xlsx")
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
inf=["VAF95%","VAF97%","VAF90%","5synergies","1synergy","4synergies"]
for i=1:length(inf)

zl010(1,:)=data.data(i,1:3);
zl010(2,:)=data.data(i,10:12);
zl010(3,:)=data.data(i,19:21);
%%%%
zl020(1,:)=data.data(i,4:6);
zl020(2,:)=data.data(i,13:15);
zl020(3,:)=data.data(i,22:24);
%%%%;
zl030(1,:)=data.data(i,7:9);
zl030(2,:)=data.data(i,16:18);
zl030(3,:)=data.data(i,25:27);

numsynl(:,1,:)=zl010;
numsynl(:,2,:)=zl020;
numsynl(:,3,:)=zl030;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zr010(1,:)=data1.data(i,1:3);
zr010(2,:)=data1.data(i,10:12);
zr010(3,:)=data1.data(i,19:21);
%%%%
zr020(1,:)=data1.data(i,4:6);
zr020(2,:)=data1.data(i,13:15);
zr020(3,:)=data1.data(i,22:24);
%%%%;
zr030(1,:)=data1.data(i,7:9);
zr030(2,:)=data1.data(i,16:18);
zr030(3,:)=data1.data(i,25:27);

numsynr(:,1,:)=zr010;
numsynr(:,2,:)=zr020;
numsynr(:,3,:)=zr030;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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
minbar=min(min(min(numsynl)))
maxbar=max(max(max(numsynl)))
ax1 = nexttile;
xx1(:,:)=numsynl(1,:,:)
b1=bar3(ax1,xx1,1,'grouped')
if minbar>10
zlim([minbar*0.98 maxbar*1.02])
else
zlim([0 maxbar])
end
set(gca,'yticklabel',yy)
title("0.8m/s")
view(-20,20);
b1(1).FaceColor='#0D4A70'
b1(2).FaceColor='#226E9C'
b1(3).FaceColor='#3C93C2'

xx2(:,:)=numsynl(2,:,:)
ax2 = nexttile;
b2=bar3(ax2,xx2,1,'grouped')
if minbar>10
zlim([minbar*0.98 maxbar*1.02])
else
zlim([0 maxbar])
end
set(gca,'yticklabel',yy)
title("1.1m/s")
view(-20,20);
b2(1).FaceColor='#0D4A70'
b2(2).FaceColor='#226E9C'
b2(3).FaceColor='#3C93C2'


xx3(:,:)=numsynl(3,:,:)
ax3 = nexttile;
b3=bar3(ax3,xx3,1,'grouped')
if minbar>10
zlim([minbar*0.98 maxbar*1.02])
else
zlim([0 maxbar])
end
set(gca,'yticklabel',yy)
title("1.45m/s")
view(-20,20);
b3(1).FaceColor='#0D4A70'
b3(2).FaceColor='#226E9C'
b3(3).FaceColor='#3C93C2'

legend({'normal', 'wide', 'wider'});
sgtitle(strcat(inf(i),"_L"),'Color','red')
set(sgtitle, 'FontSize', 20)


f2=figure()
minbar=min(min(min(numsynr)))
maxbar=max(max(max(numsynr)))
ax1 = nexttile;
xx1(:,:)=numsynr(1,:,:)
b1=bar3(ax1,xx1,1,'grouped')
set(gca,'yticklabel',yy)
if minbar>10
zlim([minbar*0.98 maxbar*1.02])
else
zlim([0 maxbar])
end
title("0.8m/s")
view(-20,20);
b1(1).FaceColor='#0D4A70'
b1(2).FaceColor='#226E9C'
b1(3).FaceColor='#3C93C2'

xx2(:,:)=numsynr(2,:,:)
ax2 = nexttile;
b2=bar3(ax2,xx2,1,'grouped')
set(gca,'yticklabel',yy)
if minbar>10
zlim([minbar*0.98 maxbar*1.02])
else
zlim([0 maxbar])
end
title("1.1m/s")
view(-20,20);
b2(1).FaceColor='#0D4A70'
b2(2).FaceColor='#226E9C'
b2(3).FaceColor='#3C93C2'

xx3(:,:)=numsynr(3,:,:)
ax3 = nexttile;
b3=bar3(ax3,xx3,1,'grouped')
set(gca,'yticklabel',yy)
if minbar>10
zlim([minbar*0.98 maxbar*1.02])
else
zlim([0 maxbar])
end
title("1.45m/s")
view(-20,20);
b3(1).FaceColor='#0D4A70'
b3(2).FaceColor='#226E9C'
b3(3).FaceColor='#3C93C2'

legend({'normal', 'wide', 'wider'});
sgtitle(strcat(inf(i),"_R"),'Color','red')
set(sgtitle, 'FontSize', 20)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(f, strcat(inf(i),"_L"),'-dpdf','-painters','-bestfit');
print(f2, strcat(inf(i),"_R"),'-dpdf','-painters','-bestfit');

end


