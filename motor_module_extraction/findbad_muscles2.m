
VAF1=95;
VAF2=97;
VAF3=90;
numsyn=5

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
for n = 1:length(source_files);
    filename=source_files(n);
    datan(1,n)=filename;
    solution=importdata(strcat(filename,".sto"));
    syndata=importdata(strcat(filename,"_19-Jul-2023.mat"));
    leftsyn=syndata.synergyOutLeft;
    Rightsyn=syndata.synergyOutRight;
    for e=1:2
        if e==1
            synergyData=Rightsyn;
       
            sided="_R";
            synData=syndata.synergyOutRight;
            
        elseif e==2
            synergyData=leftsyn;
            sided="_L";
            synData=syndata.synergyOutLeft;
        end
        %finding the number of muscle synergy needed for each leg
        temp = find([synData(:).VAF]>VAF1);
        temp2 = find([synData(:).VAF]>VAF2);
        temp3 = find([synData(:).VAF]>VAF3);
        temp4= synData(temp).VAF
        temp5= synData(1).VAF
        temp6= synData(4).VAF
        overallNSYN = temp(1);
        overallNSYN2 = temp2(1);
        overallNSYN3 = temp3(1);
        VAFof5 = temp4(1);
        VAFof1 = temp5(1);
        VAFof4 = temp6(1);
        
         datan(2,n)= overallNSYN
         datan(3,n)= overallNSYN2
        synumbers=[overallNSYN;overallNSYN2;overallNSYN3;VAFof5;VAFof1;VAFof4];
        clear temp
        clear temp2

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % making number of synergies needed excel file:
 vafname1 =  strcat(filename,'_',sided);
 vafname1=convertCharsToStrings(vafname1);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %making a file for vaf of 6 synergies:
 vafname1 =  strcat(filename,'_',sided);
 vafname1=convertCharsToStrings(vafname1);
 row1(1,2*n+e-2)=vafname1;
 ll=length(synumbers);
 row1(2:ll+1,2*n+e-2)=string(synumbers);


 %OFASym=[row,assymdata];  
%xlswrite("14530normal3.xlsx",OFASym.');
end
end
%xlswrite(strcat('badmuscles_',num2str(NoSyn),"syn"),row)
xlswrite(strcat('vaf_of_synnumbers'),row1);

