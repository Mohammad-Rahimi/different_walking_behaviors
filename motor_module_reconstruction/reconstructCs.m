tic
% list of the simulation solution files ( you also need to have extracted
% synergies from all the solution listed below, because the code compares
% them by plotting.
source_files=["08000normalsecond",
"08000wide2second",
"08000widersecond",
 "08015normal4second",
"08015widesecond",
"08015widersecond",
 "08030normal4second",
"08030wide4second",
"08030wider3second",
 "11000normalsecond",
"11000widesecond",
"11000widersecond",
 "11015normal3second",
"11015widesecond",
"11015widersecond",
 "11030normal7second",
"11030wide2second",
"11030wider5second",
"14500normalsecond",
"14500widesecond",
"14500widersecond",
"14515normalsecond",
"14515widesecond",
"14515wider2second",
"14530normalsecond",
"14530widesecond",
"14530wider2second"
];
%source_files = ["08030wider2"]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=[0
0.01
0.02
0.03
0.04
0.05
0.06
0.07
0.08
0.09
0.1
0.11
0.12
0.13
0.14
0.15
0.16
0.17
0.18
0.19
0.2
0.21
0.22
0.23
0.24
0.25
0.26
0.27
0.28
0.29
0.3
0.31
0.32
0.33
0.34
0.35
0.36
0.37
0.38
0.39
0.4
0.41
0.42
0.43
0.44
0.45
0.46
0.47
0.48
0.49
0.5
0.51
0.52
0.53
0.54
0.55
0.56
0.57
0.58
0.59
0.6
0.61
0.62
0.63
0.64
0.65
0.66
0.67
0.68
0.69
0.7
0.71
0.72
0.73
0.74
0.75
0.76
0.77
0.78
0.79
0.8
0.81
0.82
0.83
0.84
0.85
0.86
0.87
0.88
0.89
0.9
0.91
0.92
0.93
0.94
0.95
0.96
0.97
0.98
0.99
1];
theta=time*2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting up data storage
data1(2,1)="Std";
data1(3,1)="Std(xx-yy)";
data1(4,1)="mean value";
data1(5,1)="mean value reference";
data1(6,1)="Ste(xx-yy)";
data1(7,1)="sol CoA";
data1(8,1)="ref sol CoA";
data1(9,1)="R target sol";
data1(10,1)="R recont lsqr";
data1(11,1)="lsqr recons CoA";
data1(12,1)="nnlsqr recons CoA";
data1(13,1)="nnlsqr recons R";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data2(1,2)="right lsq cost";
data2(1,3)="left lsq cost";
data2(1,4)="right lsqnnM cost";
data2(1,5)="left lsqnnM cost";
data2(1,6)="right lsqnnM VAF";
data2(1,7)="left lsqnnM VAF";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% grabbing all the solutions synergy files one by one
allCrlsqnnM=[];
allCllsqnnM=[];
for n = 1:length(source_files);
    filename=source_files(n);
    data1(1,8*n-6)=strcat(filename,"_r");
    data1(1,8*n-2)=strcat(filename,"_l");
    data2(n+1,1)=filename;
    allData=importdata(strcat(filename,".sto"));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% identify which columns contain muscle ACTIVATION data
muscleActivationColumns = allData.colheaders(contains(allData.colheaders,'activation'));
muscleActivationColumns1=muscleActivationColumns;
dataMatrixAll = allData.data(:,contains(allData.colheaders,'activation'));
dataMatrixAll1=dataMatrixAll;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% removing trunk muscles from the muscle activation columns:  intobl, 
%% recabd, extobl, ercspn
remove_list=["intobl", "recabd", "extobl", "ercspn"]
for i=1:length(remove_list)
index = find(contains(muscleActivationColumns,remove_list(i)));
muscleActivationColumns(index) = [];
dataMatrixAll(:,index)=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% new muscle names based on the order we want:
muscleNames2={'piri_r','psoas_r','sart_r','tfl_r','iliacus_r','addbrev_r','addlong_r','addmagIsch_r','addmagDist_r','addmagMid_r','addmagProx_r',...
    'glmax1_r','glmax2_r','glmax3_r','glmed1_r','glmed2_r','glmed3_r',...
    'glmin1_r','glmin2_r','glmin3_r','semimem_r','bflh_r','bfsh_r','grac_r','semiten_r',...
    'vasint_r','vaslat_r','vasmed_r','recfem_r',...
    'edl_r','ehl_r','perlong_r','pertert_r','tibant_r','fdl_r','fhl_r','perbrev_r','tibpost_r','fdb_r',...
    'gaslat_r','gasmed_r','soleus_r','popli_r',...
    'piri_l','psoas_l','sart_l','tfl_l','iliacus_l','addbrev_l','addlong_l','addmagIsch_l','addmagDist_l','addmagMid_l','addmagProx_l',...
    'glmax1_l','glmax2_l','glmax3_l','glmed1_l','glmed2_l','glmed3_l',...
    'glmin1_l','glmin2_l','glmin3_l','semimem_l','bflh_l','bfsh_l','grac_l','semiten_l',...
    'vasint_l','vaslat_l','vasmed_l','recfem_l',...
    'edl_l','ehl_l','perlong_l','pertert_l','tibant_l','fdl_l','fhl_l','perbrev_l','tibpost_l','fdb_l',...
    'gaslat_l','gasmed_l','soleus_l','popli_l'};
%%%
muscleNames1=['piri_r','psoas_r','sart_r','tfl_r','iliacus_r','addbrev_r','addlong_r','addmagIsch_r','addmagDist_r','addmagMid_r','addmagProx_r',...
    'glmax1_r','glmax2_r','glmax3_r','glmed1_r','glmed2_r','glmed3_r',...
    'glmin1_r','glmin2_r','glmin3_r','semimem_r','bflh_r','bfsh_r','grac_r','semiten_r',...
    'vasint_r','vaslat_r','vasmed_r','recfem_r',...
    'edl_r','ehl_r','perlong_r','pertert_r','tibant_r','fdl_r','fhl_r','perbrev_r','tibpost_r','fdb_r',...
    'gaslat_r','gasmed_r','soleus_r','popli_r',...
    'piri_l','psoas_l','sart_l','tfl_l','iliacus_l','addbrev_l','addlong_l','addmagIsch_l','addmagDist_l','addmagMid_l','addmagProx_l',...
    'glmax1_l','glmax2_l','glmax3_l','glmed1_l','glmed2_l','glmed3_l',...
    'glmin1_l','glmin2_l','glmin3_l','semimem_l','bflh_l','bfsh_l','grac_l','semiten_l',...
    'vasint_l','vaslat_l','vasmed_l','recfem_l',...
    'edl_l','ehl_l','perlong_l','pertert_l','tibant_l','fdl_l','fhl_l','perbrev_l','tibpost_l','fdb_l',...
    'gaslat_l','gasmed_l','soleus_l','popli_l'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grab and store muscle names
temp = strfind(muscleActivationColumns,'/');
muscleNames={}
for i=1:length(temp)
    muscleNames{i} = muscleActivationColumns{i}(11:(temp{i}(3)-1));
end
% grab time data from all muscles
timeVector = allData.data(:,1);

% activation data from right leg muscles
dataMatrixR = dataMatrixAll(:,contains(muscleNames,'_r'));
% activation data from left leg muscles
dataMatrixL = dataMatrixAll(:,contains(muscleNames,'_l'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reordering muscles based on the functional order MuscleNames2

dataMatrixR2=[]
dataMatrixL2=[]
for i=1:43
dataMatrixR2(:,i) = dataMatrixR(:,find(contains(muscleNames,muscleNames2(i))));
end
for j=44:86
dataMatrixL2(:,j-43) = dataMatrixL(:,find(contains(muscleNames,muscleNames2(j)))-43);
end
dataMatrixR2=dataMatrixR2';
dataMatrixL2=dataMatrixL2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% importing reference solution synergy 
syndata1=importdata('11000normalsecond_19-Jul-2023.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reoconstructing C with different methods
%getting reference solution Ws for left and right
wr=syndata1.synergyOutRight(4).W;
wl=syndata1.synergyOutLeft(4).W;

%% least square
newCr= inv(wr'*wr)*(wr'*dataMatrixR2);
newCl= inv(wl'*wl)*(wl'*dataMatrixL2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c1=[];
c2=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% right vectorizing the problem and using matlab solver
%creating a new giant W matrix:
ll=length(newCr(1,:))*length(wr(:,1));
mm=length(newCr(1,:))*length(wr(1,:));
for ii=1:length(wr(:,1)) 
    ff=(ii-1)*length(newCr(1,:));
    bb=0;
    for hh=1:length(newCr(1,:));
        bb=hh;
    for jj=1:length(wr(1,:))
        wwr(ff+hh,bb)=wr(ii,jj);
        bb=bb+length(newCr(1,:));
    end
    end
end
%creating vectorize emg from firt time point of first muscle to last time
%point of last muscle
kk=0;
emg=[];
for j=1:length(wr(:,1))
% vectorizing emg:
kk=length(wr(:,1));
ee=dataMatrixR2(j,:);
emg=[emg,ee];
end
emgg=emg;
kk=0;
emgg=emgg';
%% solving it with matlab x = lsqnonneg(C,d)
Crlsqnonneg=lsqnonneg(wwr,emgg);
%% make c matrix again
for kkk=1:length(newCr(:,1))
    b=(kkk-1)*length(newCr(1,:))+1;
    t=(kkk)*length(newCr(1,:));
    CrlsqnnM(kkk,:)=Crlsqnonneg(b:t);
end
allCrlsqnnM=[allCrlsqnnM;CrlsqnnM];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% left vectorizing the problem and using matlab solver
%creating a new giant W matrix:
ll=length(newCl(1,:))*length(wl(:,1));
mm=length(newCl(1,:))*length(wl(1,:));
for ii=1:length(wl(:,1)) 
    ff=(ii-1)*length(newCl(1,:));
    bb=0;
    for hh=1:length(newCl(1,:));
        bb=hh;
    for jj=1:length(wl(1,:))
        wwl(ff+hh,bb)=wl(ii,jj);
        bb=bb+length(newCl(1,:));
    end
    end
end
%creating vectorize emg from firt time point of first muscle to last time
%point of last muscle
kk=0;
emgl=[];
for j=1:length(wl(:,1))
% vectorizing emg:
kk=length(wl(:,1));
ee=dataMatrixL2(j,:);
emgl=[emgl,ee];
end
emggl=emgl;
kk=0;
emggl=emggl';
% solving it with matlab x = lsqnonneg(C,d)
Cllsqnonneg=lsqnonneg(wwl,emggl);
%% (make c matrix again)
for kkk=1:length(newCr(:,1))
    b=(kkk-1)*length(newCr(1,:))+1;
    t=(kkk)*length(newCr(1,:));
    CllsqnnM(kkk,:)=Cllsqnonneg(b:t);
end
allCllsqnnM=[allCllsqnnM;CllsqnnM];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using matlab ML divide: (similar to least square)
mlCr = wr\dataMatrixR2;
mlCr = mldivide(wr,dataMatrixR2);
mlCl = wl\dataMatrixL2;
mlCl = mldivide(wl,dataMatrixL2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculating VAF for reconstructed Cs 
%% Right ( with the use of Frobenous norm)
EEr = (norm(dataMatrixR2,"fro"))^2;
eer = (norm((wr*CrlsqnnM-dataMatrixR2),"fro"))^2;
VAFlsqnnr=1-(eer/EEr);

%% left  ( with the use of Terace of (A'*A))
EEl = trace(dataMatrixL2'*dataMatrixL2);
eel = trace((wl*CllsqnnM-dataMatrixL2)'*(wl*CllsqnnM-dataMatrixL2));  
VAFlsqnnl=1-(eel/EEl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% comparing cost of different methods to least square
data2(n+1,2)=norm(wr*newCr-dataMatrixR2);
data2(n+1,3)=norm(wl*newCl-dataMatrixL2);
data2(n+1,4)=norm(wr*CrlsqnnM-dataMatrixR2);
data2(n+1,5)=norm(wl*CllsqnnM-dataMatrixL2);
% writing the VAF of reconstructed Cs into the datastorage
data2(n+1,6)=VAFlsqnnr;
data2(n+1,7)=VAFlsqnnl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% importing the target solution synergies so we can compare precvious C
%% with reference C  change the date in the file name based on the line below 
syndata=importdata(strcat(filename,"_19-Jul-2023.mat"));
synergyDataL=syndata.synergyOutLeft;    
synergyDataR=syndata.synergyOutRight;

temp1 = find([synergyDataR(:).VAF]>95);
temp2 = find([synergyDataL(:).VAF]>95);
NSYN1 = temp1(1);
NSYN2 = temp2(1);
clear temp
clear temp2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% comparing the two Cs with standard error
[r_r,p_r] = corr(syndata1.synergyOutRight(4).W, syndata.synergyOutRight(NSYN1).W);
% finding the similar Ws to compare their corresponding Cs to each other
maxr1=find(r_r(1,:)==max(r_r(1,:)));
maxr2=find(r_r(2,:)==max(r_r(2,:)));
maxr3=find(r_r(3,:)==max(r_r(3,:)));
maxr4=find(r_r(4,:)==max(r_r(4,:)));
maxr=[maxr1,maxr2,maxr3,maxr4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% right side plot
figg=figure('Renderer', 'painters', 'Position', [10 10 1600 900]);
for j=1:length(maxr)
nn=length(syndata.synergyOutRight(NSYN1).C(maxr(j),:));
xx=syndata1.synergyOutRight(4).C(j,:);
%pdgcrs=pdgCr(j,:);
lsqnnright=CrlsqnnM(j,:);
ReR=newCr(j,:);
CmlCr=mlCr(j,:);
yy=syndata.synergyOutRight(NSYN1).C(maxr(j),:);

subplot(2,2,j);
plot(time(:)/time(end)*100,xx(:),'m:','LineWidth',1);
hold on
plot(time(:)/time(end)*100,yy(:),'LineWidth',1);
hold on
plot(time(:)/time(end)*100,ReR(:),'LineWidth',1);
hold on
plot(time(:)/time(end)*100,CmlCr(:),'LineWidth',1);
hold on
plot(time(:)/time(end)*100,lsqnnright(:),'LineWidth',1);
    xlabel('Stride (%)','FontWeight','b');
    ylabel('activation','FontWeight','b');
    legend(' C of reference sol',...
        'c of target solution(title)','least square',...
        'matlab least square',"matlab nonnegative least square");
    %ylim([-30 30])
    title([strcat(filename,'R Cs') num2str(j) num2str(maxr(j))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculating the center of activity
xxCoa=(atan2((xx*sin(theta)),(xx*cos(theta))));
yyCoa=(atan2((yy*sin(theta)),(yy*cos(theta))));
ReRCoa=(atan2((ReR*sin(theta)),(ReR*cos(theta))));
RennRCoa=(atan2((lsqnnright*sin(theta)),(lsqnnright*cos(theta))));%%lsqnnright
if xxCoa<0
    xxCoa=xxCoa+2*pi;
end
if yyCoa<0
    yyCoa=yyCoa+2*pi;
end
if ReRCoa<0
    ReRCoa=ReRCoa+2*pi;
end
if RennRCoa<0
    RennRCoa=RennRCoa+2*pi;
end
xxCoa=xxCoa/(2*pi);
yyCoa=yyCoa/(2*pi);
ReRCoa=ReRCoa/(2*pi);
RennRCoa=RennRCoa/(2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% right standard error and correlation between Cs
% calculating standard error between two C and standar deviation of each Cs
% because If time series xx is the similar to time series yy then the
% variance of xx-yy should be less than the variance of xx.
sterror1=sqrt(((xx-yy)*(xx-yy)')/(nn-1));
StE1=std(xx-yy);
sd1= std(yy);
mr1=mean(xx);
mrr1= mean(yy);
[CRr,CPr] = corr(xx',yy');
[CRER,CPER]=corr(xx',ReR');
[CRERnn,CPERnn]=corr(xx',lsqnnright');%lsqnnright
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% writing everything into a excel file

data1(2,8*n-7+j)=sd1;
data1(3,8*n-7+j)=StE1;
data1(4,8*n-7+j)=mrr1;
data1(5,8*n-7+j)=mr1;
data1(6,8*n-7+j)=sterror1;
data1(7,8*n-7+j)=yyCoa;
data1(8,8*n-7+j)=xxCoa;
data1(9,8*n-7+j)=CRr;
data1(10,8*n-7+j)=CRER;
data1(11,8*n-7+j)=ReRCoa;
data1(12,8*n-7+j)=RennRCoa;
data1(13,8*n-7+j)=CRERnn;

end
%exportgraphics(figg,strcat("000000",filename,"ReCons_Cs_R",'.jpg'),'Resolution',1000)




%xlswrite(strcat("14500normal3Syn",filename,"nSyn",num2str(NSYN1),"_R"),[r_r;aa;p_r])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finding similar Ws based on their R value and compare the corresponding Cs
[r_l,p_l] = corr(syndata1.synergyOutLeft(4).W, syndata.synergyOutLeft(NSYN2).W);
bb=zeros(NSYN2);
% finding the similar Ws to compare their corresponding Cs to each other
maxl1=find(r_l(1,:)==max(r_l(1,:)));
maxl2=find(r_l(2,:)==max(r_l(2,:)));
maxl3=find(r_l(3,:)==max(r_l(3,:)));
maxl4=find(r_l(4,:)==max(r_l(4,:)));
maxl=[maxl1,maxl2,maxl3,maxl4];

figg=figure('Renderer', 'painters', 'Position', [10 10 1600 900]);
for j=1:length(maxl)
nn=length(syndata.synergyOutLeft(NSYN2).C(maxl(j),:));
xx=syndata1.synergyOutLeft(4).C(j,:);
ReL=newCl(j,:);
%pdgcls=pdgCl(j,:);
CmlCl=mlCl(j,:);
lsqnnleft=CllsqnnM(j,:);
yy=syndata.synergyOutLeft(NSYN2).C(maxl(j),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot left
subplot(2,2,j);
plot(time(:)/time(end)*100,xx(:),'m:','LineWidth',1);
hold on
plot(time(:)/time(end)*100,yy(:),'LineWidth',1);
hold on
plot(time(:)/time(end)*100,ReL(:),'LineWidth',1);
hold on
plot(time(:)/time(end)*100,CmlCl(:),'LineWidth',1);
hold on
plot(time(:)/time(end)*100,lsqnnleft(:),'LineWidth',1);
    xlabel('Stride (%)','FontWeight','b');
    ylabel('activation','FontWeight','b');
    legend(' C of reference sol',...
        'c of target solution(title)','least square sol',...
        'matlab least square sol',"nonnegative least square sol");
    %ylim([-30 30])
    title([strcat(filename,' L Cs ') num2str(j) num2str(maxl(j))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculating center of the activity
xxCoa=(atan2((xx*sin(theta)),(xx*cos(theta))));
yyCoa=(atan2((yy*sin(theta)),(yy*cos(theta))));
ReLCoa=(atan2((ReL*sin(theta)),(ReL*cos(theta))));
RennLCoa=(atan2((lsqnnleft*sin(theta)),(lsqnnleft*cos(theta))));%lsqnnleft
if xxCoa<0
    xxCoa=xxCoa+2*pi;
end
if yyCoa<0
    yyCoa=yyCoa+2*pi;
end
if ReLCoa<0
    ReLCoa=ReLCoa+2*pi;
end
if RennLCoa<0
    RennLCoa=RennLCoa+2*pi;
end
xxCoa=xxCoa/(2*pi);
yyCoa=yyCoa/(2*pi);
ReLCoa=ReLCoa/(2*pi);
RennLCoa=RennLCoa/(2*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Left standard error and correlation between two Cs
% calculating standard error between two C and standar deviation of each Cs
% because If time series xx is the similar to time series yy then the variance
% of xx-yy should be less than the variance of xx.
sterror1=sqrt(((xx-yy)*(xx-yy)')/(nn-1));
StE1=std(xx-yy);
sd1= std(yy);
mr1=mean(xx);
mrr1= mean(yy);
[CRl,CPl] = corr(xx',yy');
[CREL,CPEL] = corr(xx',ReL');
[CRELnn,CPELnn]=corr(xx',lsqnnleft');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% writing everything into a excel file
data1(2,8*n-3+j)=sd1;
data1(3,8*n-3+j)=StE1;
data1(4,8*n-3+j)=mrr1;
data1(5,8*n-3+j)=mr1;
data1(6,8*n-3+j)=sterror1;
data1(7,8*n-3+j)=yyCoa;
data1(8,8*n-3+j)=xxCoa;
data1(9,8*n-3+j)=CRl;
data1(10,8*n-3+j)=CREL;
data1(11,8*n-3+j)=ReLCoa;
data1(12,8*n-3+j)=RennLCoa;
data1(13,8*n-3+j)=CRELnn;
end
%exportgraphics(figg,strcat(filename,"ReCons_Cs_L",'.jpg'),'Resolution',1000)
%xlswrite(strcat("14500normal3Syn",filename,"nSyn",num2str(NSYN1),"_L"),[r_l;bb;p_l])
end
xlswrite("compareCs_reconstfinalJan",data1)
xlswrite("compareCs_costsfinalJan",data2)
xlswrite("all_recons_C_r_finalJan",allCrlsqnnM)
xlswrite("all_recons_C_l_finalJan",allCllsqnnM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc

