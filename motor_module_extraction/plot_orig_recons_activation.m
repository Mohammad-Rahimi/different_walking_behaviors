%newfilename='test'
NoSyn=5
filename='14500normal'
solution=importdata(strcat(filename,".sto"))
syndata=importdata(strcat(filename,"_30-Jan-2023.mat"));
timeVector=solution.data(:,1)
muscleActivationColumns = solution.colheaders(contains(solution.colheaders,'activation'));
%%% removing trunk muscles from the muscle activation columns:  intobl, recabd, extobl, ercspn
remove_list=["intobl", "recabd", "extobl", "ercspn"]
for i=1:length(remove_list)
index = find(contains(muscleActivationColumns,remove_list(i)))
muscleActivationColumns(index) = []
end

% grab and store muscle names
temp = strfind(muscleActivationColumns,'/');
for i=1:length(temp)
    muscleNames{i} = muscleActivationColumns{i}(11:(temp{i}(3)-1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
leftsynergy=syndata.synergyOutLeft
rightsynergy=syndata.synergyOutRight
%temp = find([synergyData(:).VAF]>95);
%temp2 = find([synergyData(:).VAF]>97);
% NSYN1 = temp(1);
% NSYN2 = 6;
% clear temp
% clear temp2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leftsyn=syndata.synergyOutLeft(NSYN1)
%Rightsyn=syndata.synergyOutRight(NSYN1)

for k=1:2
    if k==1
      synergydata=leftsynergy
      temp = find([synergydata(:).VAF]>95);
      NSYN1 = temp(1);
      leftsyn=synergydata(NSYN1)
      %NSYN2 = 6;
      clear temp
      %clear temp2
      

      ws=leftsyn.W
      cs=leftsyn.C
      newfilename='_L'
      synergyData=leftsyn
      names = muscleNames(44:end);
    else
      synergydata=rightsynergy
      temp = find([synergydata(:).VAF]>95);
      NSYN1 = temp(1);
      rightsyn=synergydata(NSYN1)
      %NSYN2 = 6;
      clear temp
      ws=rightsyn.W
      cs=rightsyn.C
      newfilename='_R'
      synergyData=rightsyn
      names = muscleNames(1:43);
    end
    i=1;
for m=1:length(names)
        temp1=[];
        temp1=synergyData.VAFmus(m);
        if temp1<75
            xx(i,:)=synergyData.emgOrig(m,:);
            yy(i,:)=synergyData.emgRecon(m,:);
            o(i)=names(m);
            i=i+1;
        end
        totalnumber=i-1
        

end
time=timeVector;
figg=figure('Renderer', 'painters', 'Position', [10 10 1600 900]);
for j=1:totalnumber;
newfname=erase(filename,".sto")
 subplot(floor((totalnumber/2))+1,2,j)
  
    
        
            
           
            hold on ;
            box on;
            plot(time(:)/time(end)*100,xx(j,:),'LineWidth',1);
            plot(time(:)/time(end)*100,yy(j,:),'LineWidth',1);
            xlabel('Stride (%)','FontWeight','b');
            ylabel('activation' ,'FontWeight','b');
            title(string(o(j)))
            legend(["origin"  "recons"])
end
end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
    