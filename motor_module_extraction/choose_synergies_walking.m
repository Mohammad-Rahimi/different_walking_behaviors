function [nsyn] = choose_synergies_walking(fname,synergyData,CI,chooseFlag,emgNames,dataMatrix,timeVector)


[overallNSYN,overallNSYN2,ciNSYN,musNSYN] = make_datafile(fname,synergyData,CI,emgNames);

nsyn1 = ciNSYN;
%nsyn2 = max([overallNSYN,musNSYN]);
nsyn2 = overallNSYN;
nsyn3 = overallNSYN2;
% choose number of synergies
if chooseFlag == 1 % use overall VAF with CI only
    nsyn = nsyn1;
else % use max of all VAFs (minus background condition)
    nsyn = nsyn2;
end


makeVAFplots(fname,synergyData,CI,emgNames,nsyn1,nsyn2,nsyn3,dataMatrix,timeVector);


function [overallNSYN,overallNSYN2,ciNSYN,musNSYN,condNSYN] = make_datafile(fname,synergyData,CI,emgNames)
ciNSYN=10;

% Make file with numerical data
% nSynTotal = size(emgNames,1);
nSynTotal = size(synergyData,2);

% open file
fid = fopen([fname '_' date '.dat'],'w');

% create headers
fprintf(fid,'nSyn\t');

for i=1:nSynTotal
    fprintf(fid,'%i\t',i);
end
fprintf(fid,'\n');

%*******************%
%--- overall VAF ---%
%*******************%
fprintf(fid,'Overall VAF\t');
for i=1:nSynTotal
    fprintf(fid,'%5.2f\t',synergyData(i).VAF);
end
fprintf(fid,'\n');

temp = find([synergyData(:).VAF]>95);
temp2 = find([synergyData(:).VAF]>97);
overallNSYN = temp(1);
overallNSYN2 = temp2(1);
clear temp
clear temp2
%**********************************%
%--- overall VAF CI lower bound ---%
%**********************************%
% fprintf(fid,'Overall VAF CI\t');
% 
% DATA = synergyData;
% 
% % how many replicates?
% if ~isfield(DATA,'sortedVAF_original')
%     nreps = 100;
% else
%     nreps = size(DATA(1).sortedVAFmus_original,2);
% end
% 
% temp = (100 - CI)/(2*100);
% up = ceil((1 - temp)*nreps);
% bot = ceil(temp*nreps);
% 
% for i=1:nSynTotal
%     if ~isfield(DATA,'sortedVAF_original')
%         lb(i) = DATA(i).VAF;
%     else
%     lb(i) = DATA(i).sortedVAF_original(bot);
%     end
%     fprintf(fid,'%5.2f\t',lb(i));
% end
% fprintf(fid,'\n');
% 
% temp = find(lb>90);
% if ~isempty(temp)
%     ciNSYN = temp(1);
% else
%     ciNSYN = NaN;
% end
% clear temp

%******************%
%--- muscle VAF ---%
%******************%
musNSYN = [];
for m=1:length(emgNames)
    fprintf(fid,'%s\t',emgNames{m});
    for i=1:nSynTotal
        fprintf(fid,'%5.2f\t',synergyData(i).VAFmus(m));
        if isempty(musNSYN)
            temp = find([synergyData(i).VAFmus]<75);
            if isempty(temp)
                musNSYN = i;
            end
            clear temp
        end
    end
    fprintf(fid,'\n');
end


fclose(fid);

function makeVAFplots(fname,synergyData,CI,emgNames,nsyn1,nsyn2,nsyn3,dataMatrix,timeVector)

% Plot VAF
f=figure('units','inches','position',[0 -1 11 8.5]);


% COLORS FOR PLOTTING
%                  1          2          3          4          5          6          7        8           9        10            11        12
%                red     fade red     green    fade green    blue     fade blue   black     gray       yellow    orange        plum      lime
color_palette= [255 0 0; 255 150 150; 0 114 54; 153 255 102; 0 0 255; 44 176 207; 0 0 0; 155 155 155; 255 240 0; 255 185 0; 153 51 102; 171 218 77]*(1/255);
color=color_palette([1 10 4 6 11 3 5 2 12 9 8 1 10 4 6 11 3 2 5 12 1 10 4 6 11 3 5 2 12 9 8 1 10 4 6 11 3 2 5 12 1 10 4 6 11 3 2 5 12],:);


%***************************************%
% OVERALL VAF WITH CONFIDENCE INTERVALS %
%***************************************%

DATA = synergyData;
% how many replicates?
if ~isfield(DATA,'sortedVAF_original')
    nreps = 100;
else
    nreps = size(DATA(1).sortedVAFmus_original,2);
end
temp = (100 - CI)/(2*100);
up = ceil((1 - temp)*nreps);
bot = floor(temp*nreps);
DATA = synergyData(:);

% For overall VAF
for i=1:length(DATA)
    
    if ~isfield(DATA,'sortedVAF_original')
        lb(i) = DATA(i).VAF;
        ub(i) = DATA(i).VAF; 
        lb2(i) = 0;
        ub2(i) = 0;
        x2(i) = 0;
    else
        lb(i) = DATA(i).sortedVAF_original(bot);
        ub(i) = DATA(i).sortedVAF_original(up);
        lb2(i) = DATA(i).sortedVAF_shuffled(bot);
        ub2(i) = DATA(i).sortedVAF_shuffled(up);
        x2(i) = (ub2(i)-lb2(i))/2 + lb2(i);
    end
    
    
end

subplot(2,2,1)
hold on
errorbar(1:length(DATA),[DATA(:).VAF], ([DATA(:).VAF]-lb), (ub-[DATA(:).VAF]));
%errorbar(1:length(DATA), x2, (x2-lb2), (ub2-x2),'color','r');
xx=[0 10];
yy= [95 95];
line(xx,yy,'Color','red','LineStyle','--');
xlim([0 size(synergyData,2)]);
title([strcat(fname,' Overal VAF with') num2str(CI) '% CI']);
legend('real','overal VAF 95%','location','best');


%************%
% MUSCLE VAF %
%************%

temp = [synergyData(:).VAFmus];
[nmus,nsyn] = size(temp);
subplot(2,2,2);
for m=1:nmus
    temp1=[];
    for k=1:nsyn
        temp1(k)=synergyData(k).VAFmus(m);
        axis([0 nsyn 0 100]);
        set(gca,'Xtick',1:nmus);
        xtickangle(90)
    end
    %plot([0:nsyn],[0,temp1],'LineWidth',2,'Color',color(m,:))
    plot([0:nsyn],[0,temp1],'LineWidth',2);
    hold on;
end
title(['VAF each muscle']);
set(gca,'Ytick',[25 50 75 100]);
%legend(emgNames);
xlabel('number of components');
ylabel('VAF (%)');

%***************%
% CONDITION VAF %
%***************%

% placeholder for gait cycle bin VAF

% index=1;
% for i=1:size(synergyData,2)
%     legend_names{i} = [num2str(i)];
% end
% 
% VAFcond=[];
% for b=1:numbins %num of bins
%     for k=1:size(synergyData,2)
%         
%         VAFcond(:,:,:,k)=reshape(synergyData(k).VAFcond,numtrials,numconds,numbins);
%         
%         subplot(2,4,index+4)
%         hold on
%         plot(dirs,nanmean(VAFcond(:,:,b,k)),'LineWidth',2,'Color',color(k,:))
%         axis([0 330 0 100])
%         set(gca,'Ytick',[25 50 75 80 85 90 95 100])
%         set(gca,'Xtick',[0 90 180 270])
%         
%     end
%     
%     index=1+index;
%     
%     if b==1
%         title('BKGD')
%         legend(legend_names(1:10),'Location','SW');
%     elseif b==2
%         title('APR1')
%     elseif b==3
%         title('APR2')
%     elseif b==4
%         title('APR3')
%     end
%     
%     
% end


subplot(2,2,3);
text(1.2,1.05,[fname ' (' date ')'],'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',18,'Interpreter','none')
% text(1.2,-0.25,['nsyn = ' num2str(nsyn1) ', VAF CI only'],'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12,'Interpreter','none');
text(1.2,-0.15,['nsyn = ' num2str(nsyn2)],'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',12,'Interpreter','none');
 X=[1:2];
 name=["VAF95%";"VAF97%"];
    Y=[nsyn2,nsyn3];
    bar(X,Y,'green');
    ylabel('number of synergies','FontWeight','b');
    set(gca,'xticklabel',name);
     xtickangle(90)
    %title(newfilename)
% Save current figure
orient tall;
print(f,'-depsc',[fname '_' date '.eps']);
newfname=erase(fname,".sto")
exportgraphics(f,strcat(newfname,'.jpg'),'Resolution',500)
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%
figg=figure('Renderer', 'painters', 'Position', [10 10 1600 900]);
    subplot(2,2,1)
    %hold on; box on;
    name = emgNames(1:nmus);
    xx=[];
    yy=[];
    i=1;
    for m=1:nmus
        temp1=[];
        
        temp1=synergyData(nsyn2).VAFmus(m);
        if temp1<75
            xx(i)=temp1;
            o(i)=emgNames(m);
            i=i+1;
                end
    end
    disp("this is the muscles below VAF 75% for overal VAF 95%")
    VAF95=o'
    %plot([0:nsyn],[0,temp1],'LineWidth',2,'Color',color(m,:))
    X=[1:length(xx)];
    bar(X,xx,'yellow');
    ylabel('VAF (%)','FontWeight','b');
    set(gca,'xticklabel',o);
     xtickangle(90)
    title(fname," muscles below 75% for overal VAF 95%");
    
    subplot(2,2,2);
    xx=[];
    o={};
    i=1;
    for m=1:nmus
        temp1=[];
        
        temp1=synergyData(nsyn3).VAFmus(m);
        if temp1<75
            xx(i)=temp1;
            o(i)=emgNames(m);
            i=i+1;
        end
    end
    disp("this is the muscles below VAF 75% for overal VAF 97%")
    VAF97=o'
    %plot([0:nsyn],[0,temp1],'LineWidth',2,'Color',color(m,:))
    X=[1:length(xx)];
    bar(X,xx,'yellow');
    ylabel('VAF (%)','FontWeight','b');
    set(gca,'xticklabel',o);
     xtickangle(90)
    title(fname," muscles below 75% for overal VAF 97%");
    subplot(2,2,3);
    o={};
    time=timeVector;
    i=1
    for m=1:nmus
        temp1=[];
        
        temp1=synergyData(nsyn2).VAFmus(m);
        if temp1<75
            xx=dataMatrix(m,:);
            o(i)=emgNames(m);
            hold on ;
            box on;
            plot(time(:)/time(end)*100,xx,'LineWidth',1);
            xlabel('Stride (%)','FontWeight','b');
            ylabel('activation' ,'FontWeight','b');
            
            ylim([0 0.5]);
            title(fname,'activation of muscles below VAF% with overal VAF 95% ');
            i=i+1;
        end
        legend(o,'Location','NorthEastOutside');
    end
subplot(2,2,4)
    o={};
    xx=[];
    time=timeVector;
    i=1;
    for m=1:nmus
        temp1=[];
        
        temp1=synergyData(nsyn3).VAFmus(m);
        if temp1<75
            xx=dataMatrix(m,:);
            %if max(xx)<0.1
            o(i)=emgNames(m);
            hold on ;
            box on;
            plot(time(:)/time(end)*100,xx,'LineWidth',1);
            xlabel('Stride (%)','FontWeight','b');
            ylabel('activation' ,'FontWeight','b');
            
            ylim([0 0.5]);
            title(fname,'activation of muscles below VAF 75% with overal VAF 97%')
            i=i+1;
            %end
        end
        legend(o,'Location','NorthEastOutside');
        
    end
newfname=erase(fname,".sto")
exportgraphics(figg,strcat(newfname,"_activation",'.jpg'),'Resolution',500)
disp("this is the muscles below VAF 75% for overal VAF 95% whose max activation is less than 10%")
vaf9510=o'
i=1;
for m=1:nmus
        temp3=[];
        temp3=synergyData(nsyn3).VAFmus(m);
        if temp3<75
            xx=dataMatrix(m,:);
            if max(xx)<0.1
            q(i)=emgNames(m);
            i=i+1;
            
            end
        end
end
disp("this is the muscles below VAF 75% for overal VAF 97% whose max activation is less than 10%")
 vaf9710=q'
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % writing muscle below 75% in an excel file : VAF97 vaf9710 vaf9510
 listt=["VAF97", "vaf9710","vaf95", "vaf9510",];
 listtt={string(VAF97), string(vaf9710),string(VAF95), string(vaf9510)};
 %if contains(VAF97,'_r');
 %    c="_r";
 %else
 %    c="_l";
 %end
 ffname=erase(fname,".sto")
 for i=1:4   
 termNamestr =  strcat(ffname,'_',listt(i));
 termNamestr=convertCharsToStrings(termNamestr);
 row(1,i)=termNamestr;
 ll=length(listtt{i});
 row(2:ll+1,i)=string(listtt{i});
 %OFASym=[row,assymdata];  
%xlswrite("14530normal3.xlsx",OFASym.');
 end
xlswrite(strcat(ffname),row)
row=[]
