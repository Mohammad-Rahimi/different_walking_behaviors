
clear all
close all

%fileLocation = 'C:\Users\Mohammad\Desktop\New folder';
filename = '11000normalsecondsecond.sto';
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
dataMatrixR2=dataMatrixR2';
dataMatrixL2=dataMatrixL2';
% % normalize to max across conditions % I commented this section out
% because these are simulation data and all simulations are on the same
% scale (i.e., muscle activity can vary from 0 to 1). We usually have to do
% this for EMG data to put them on a "fake" 0 to 1 scale. 
% emgMaxValueR = max(dataMatrixR,[],2);
% dataMatrixR = dataMatrixR./(emgMaxValueR*ones(1,size(dataMatrixR,2)));% 
% emgMaxValueL = max(dataMatrixL,[],2);
% dataMatrixL = dataMatrixL./(emgMaxValueL*ones(1,size(dataMatrixL,2)));

%dataMatrixR2=dataMatrixR2';
% extract synergies
disp('extracting synergies full set')

[synergyOutRight] = synergyExtraction_Jessica(dataMatrixR2,10,1);
[synergyOutLeft] = synergyExtraction_Jessica(dataMatrixL2,10,1);

% disp('bootstrapping synergies full set');
% nboot = 100;
% [synergyOutRight] = bootstrap_VAF_Jessica(synergyOutRight, 8, nboot);
% [synergyOutLeft] = bootstrap_VAF_Jessica(synergyOutLeft, 8, nboot);

CI = 95;
chooseFlag = 2; % choose number of synergies based on 1:CI, 2:max of all
% inputs for this: fname,synergyData,CI,chooseFlag,emgNames
[nsynR_postureFull] = choose_synergies_walking([strcat(filename,'_r')],synergyOutRight,CI,chooseFlag,muscleNames2(1:43),dataMatrixR2,timeVector);
[nsynL_postureFull] = choose_synergies_walking([strcat(filename,'_l')],synergyOutLeft,CI,chooseFlag,muscleNames2(44:end),dataMatrixL2,timeVector);

emgNames=[];
newfilename=erase(filename,".sto")
save([newfilename '_' date],'dataMatrixR','dataMatrixL','synergyOutRight','synergyOutLeft','emgNames');
%% compare synergies

% need to talk to Hannah about what code we are currently using to compare
% synergies from different conditions, subjects, etc. 


