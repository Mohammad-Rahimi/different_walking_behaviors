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
data2(1,1)="Max w Rs";

data3(1,1)="max C Rs";
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


for n = 1:length(source_files);
    filename=source_files(n);
    syndata1=importdata("11000normalsecond_19-Jul-2023.mat");
    syndata2=importdata(strcat(filename,"_19-Jul-2023.mat"));
   
    data2(n+1,1)=filename;
    data3(n+1,1)=filename;
    synergyData1=syndata1.synergyOutLeft;  
    %synergyData1=syndata1.synergyOutLeft;
    synergyData2=syndata2.synergyOutLeft;

temp1 = find([synergyData1(:).VAF]>95);
temp2 = find([synergyData2(:).VAF]>95);
NSYN1 = temp1(1);
NSYN2 = temp2(1);
clear temp
clear temp2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r_r,p_r] = corr(synergyData1(NSYN1).W, synergyData2(NSYN2).W);
%[r_rc,p_rc] = corr(synergyData1(NSYN1).C, synergyData2(NSYN2).C)

for i=1:NSYN1
    maxRs(i)=max(r_r(i,:));
    maxr1=find(r_r(1,:)==max(r_r(1,:)));
    maxr2=find(r_r(2,:)==max(r_r(2,:)));
    maxr3=find(r_r(3,:)==max(r_r(3,:)));
    maxr4=find(r_r(4,:)==max(r_r(4,:)));
    maxr=[maxr1,maxr2,maxr3,maxr4];
    data2(n+1,i+1)=max(r_r(i,:));
    yy=synergyData2(NSYN2).C(maxr(i),:);
    xx= synergyData1(NSYN1).C(i,:);
    data3(1,i+1)="max"
    data3(n+1,i+1)=max(synergyData2(NSYN2).C(maxr(i),:)');
    data3(1,i+5)="R"
    data3(n+1,i+5)=corr(synergyData1(NSYN1).C(i,:)',synergyData2(NSYN2).C(maxr(i),:)');
    data3(1,i+9)="CoA"
    yyCoa=(atan2((yy*sin(theta)),(yy*cos(theta))));
    if yyCoa<0
       yyCoa=yyCoa+2*pi;
    end
    yyCoa=yyCoa/(2*pi);
    data3(n+1,i+9)=yyCoa;
    data3(1,i+13)="dwt";
    data3(n+1,i+13)=dtw(xx',yy');

end

end
xlswrite("w_Left_comparison_to_110test",data2)
xlswrite("C_Left_comparison_to_110test",data3)