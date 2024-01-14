
clear all
close all
excelL(1,1)="left"
excelR(1,1)="Right"
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
%fileLocation = 'C:\Users\Mohammad\Desktop\New folder';
filename = '11000normalsecond.sto';
%outputName = 'ttest'; 
muscles = []; % if you want to extract synergies from a reduced set of muscles, you can list the columns here

% load in .sto file from Moco 
%checked th last header is 21
allData = importdata(filename); %confirm that the last header in your moco output file is at line 23, otherwise change

% identify which columns contain muscle ACTIVATION data
muscleActivationColumns = allData.colheaders(contains(allData.colheaders,'activation'));
muscleActivationColumns1=muscleActivationColumns;
dataMatrixAll = allData.data(:,contains(allData.colheaders,'activation'));
dataMatrixAll1=dataMatrixAll;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% removing trunk muscles from the muscle activation columns:  intobl, recabd, extobl, ercspn
remove_list=["intobl", "recabd", "extobl", "ercspn"]
for i=1:length(remove_list)
index = find(contains(muscleActivationColumns,remove_list(i)));
muscleActivationColumns(index) = [];
dataMatrixAll(:,index)=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% new muscle names based on the order we want:

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grab and store muscle names
temp = strfind(muscleActivationColumns,'/');
for i=1:length(temp)
    muscleNames{i} = muscleActivationColumns{i}(11:(temp{i}(3)-1));
end

% grab ACTIVATION data from all muscles
timeVector = allData.data(:,1);

% activation data from right leg muscles
dataMatrixR = dataMatrixAll(:,contains(muscleNames,'_r'));
% activation data from left leg muscles
dataMatrixL = dataMatrixAll(:,contains(muscleNames,'_l'));
dataMatrixLL=dataMatrixL';
%dataMatrixR=(dataMatrixR)';
%dataMatrixL=(dataMatrixL)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reordering muscles based on the functional order MuscleNames2
for i=1:43
dataMatrixR2(:,i) = dataMatrixR(:,find(contains(muscleNames,muscleNames2(i))));
end
for j=44:86
dataMatrixL2(:,j-43) = dataMatrixL(:,find(contains(muscleNames,muscleNames2(j)))-43);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(source_files)
    muscles = []; % if you want to extract synergies from a reduced set of muscles, you can list the columns here

% load in .sto file from Moco 
%checked th last header is 21
allData = importdata(strcat(source_files(ii),".sto")); %confirm that the last header in your moco output file is at line 23, otherwise change
excelL(1,ii+1)=source_files(ii)
excelR(1,ii+1)=source_files(ii)
% identify which columns contain muscle ACTIVATION data
muscleActivationColumns = allData.colheaders(contains(allData.colheaders,'activation'));
muscleActivationColumns1=muscleActivationColumns;
dataMatrixAll = allData.data(:,contains(allData.colheaders,'activation'));
dataMatrixAll1=dataMatrixAll;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% removing trunk muscles from the muscle activation columns:  intobl, recabd, extobl, ercspn
remove_list=["intobl", "recabd", "extobl", "ercspn"]
for i=1:length(remove_list)
index = find(contains(muscleActivationColumns,remove_list(i)));
muscleActivationColumns(index) = [];
dataMatrixAll(:,index)=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% new muscle names based on the order we want:

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grab and store muscle names
temp = strfind(muscleActivationColumns,'/');
for i=1:length(temp)
    muscleNames{i} = muscleActivationColumns{i}(11:(temp{i}(3)-1));
end

% grab ACTIVATION data from all muscles
timeVector = allData.data(:,1);

% activation data from right leg muscles
dataMatrixR = dataMatrixAll(:,contains(muscleNames,'_r'));
% activation data from left leg muscles
dataMatrixL = dataMatrixAll(:,contains(muscleNames,'_l'));
dataMatrixLL=dataMatrixL';
%dataMatrixR=(dataMatrixR)';
%dataMatrixL=(dataMatrixL)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reordering muscles based on the functional order MuscleNames2
for i=1:43
dataMatrixR22(:,i) = dataMatrixR(:,find(contains(muscleNames,muscleNames2(i))));
end
for j=44:86
dataMatrixL22(:,j-43) = dataMatrixL(:,find(contains(muscleNames,muscleNames2(j)))-43);
end
for jj=1:length(dataMatrixL22(1,:))
    excelL(jj+1,ii+1)=corr(dataMatrixL2(:,jj),dataMatrixL22(:,jj))
    excelR(jj+1,ii+1)=corr(dataMatrixR2(:,jj),dataMatrixR22(:,jj))
end
end
xlswrite("R_activations_L",excelL)  
xlswrite("R_activations_R",excelR)  
    
    

