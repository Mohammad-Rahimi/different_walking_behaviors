newfilename='test'
source_files=[%"08000normalsecond",
% "08000wide2second",
% "08000widersecond",
%  "08015normal4second",
% "08015widesecond",
% "08015widersecond",
%  "08030normal4second",
% "08030wide4second",
% "08030wider3second",
  "11000normalsecond",
% "11000widesecond",
% "11000widersecond",
%  "11015normal3second",
 "11015widesecond"
% "11015widersecond",
%  "11030normal7second",
% "11030wide2second",
% "11030wider5second",
% "14500normalsecond",
% "14500widesecond",
% "14500widersecond",
% "14515normalsecond",
% "14515widesecond",
% "14515wider2second",
% "14530normalsecond",
% "14530widesecond",
% "14530wider2second"
];
%source_files=["14500normal"]
for n = 1:length(source_files);
    filename=source_files(n);
    solution=importdata(strcat(filename,".sto"));
    syndata=importdata(strcat(filename,"_19-Jul-2023.mat"));
    %NoSyn=5
    timeVector=solution.data(:,1)
    muscleActivationColumns = solution.colheaders(contains(solution.colheaders,'activation'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    muscleNames={'piri_r','psoas_r','sart_r','tfl_r','iliacus_r','addbrev_r','addlong_r','addmagIsch_r','addmagDist_r','addmagMid_r','addmagProx_r',...
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     leftsyn=syndata.synergyOutLeft
%     Rightsyn=syndata.synergyOutRight
    %temp = find([synergyData(:).VAF]>95);
    %temp2 = find([synergyData(:).VAF]>97);
    %NSYN1 = temp(1);
    %NSYN2 = 6;
    clear temp
    clear temp2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     leftsyn=syndata.synergyOutLeft(NoSyn)
%     Rightsyn=syndata.synergyOutRight(NoSyn)
            synergyData=syndata.synergyOutLeft
            
            %ws=Rightsyn.W
            %cs=Rightsyn.C
            newfilename='_L'
            name = erase(muscleNames(1:43),"_L");
            %synergyData=Rightsyn
        %origi=leftsyn.emgOrig;
        time=timeVector;
        temp = find([synergyData(:).VAF]>95);
        NoSyn = temp(1);
        %NoSyn=NSYN1
        %NSYN2 = 6;
        clear temp
        %clear temp2

            %ws=leftsyn.W
            %cs=leftsyn.C
            syn=synergyData(NoSyn)
            ws=syn.W;
            %wss(:,:,n)=syn.W
            cs=syn.C;
            figg=figure('Renderer', 'painters', 'Position', [10 10 1600 900]);
            newfname=erase(filename,".sto")
            %title(strcat(newfname,newfilename,' Nsyn='),NoSyn);
            for i= 1:NoSyn
                subplot(NoSyn,2,2*i-1)
                %hold on; box on;

                X=[1:43];
                Y=ws(:,i);
                bar(X,Y,'FaceColor',"#FF4500")
                ylabel('Activation','FontWeight','b')
                set(gca,'xtick',1:length(name))
                if i==NoSyn
                set(gca,'xtick',1:length(name))
                set(gca,'xticklabel',name(:),'FontSize',7)
                xtickangle(90)
                
                end
                if i==1
                title(strcat("Ws composition",newfname,newfilename));
                end
                
                %newfilename=(erase(targetsolution,".sto"));
                %saveas(figg,strcat(newfilename,'.jpg'))



                subplot(NoSyn,2,2*i)

                plot(time(:)/time(end)*100,cs(i,:),'-',"color","#0000CD",'LineWidth',1)
                xlabel('Stride (%)','FontWeight','b')
                ylabel('activation','FontWeight','b')
                %legend('cc', 'Avg step width','l-double support','r-double support')
                %ylim([-30 30])
                if i==1
                title(strcat(newfname,newfilename,' Nsyn= ',num2str(NoSyn)));
                end
            end
            newfname=erase(filename,".sto");
                orient(figg,'landscape')
                exportgraphics(figg,strcat(newfname,newfilename,"nSyn",num2str(NoSyn),"VAF95",'.jpg'),'Resolution',500)
                exportgraphics(figg,strcat(newfname,newfilename,"nSyn",num2str(NoSyn),"VAF95",'.pdf'),'Resolution',1600)
                print(figg, strcat(newfname,newfilename,"nSyn1"),'-dpdf','-painters','-bestfit');
                [r,p] = corr(synergyData(NoSyn).W, synergyData(NoSyn).W)
                aa=zeros(NoSyn);
                xlswrite(strcat(newfname,newfilename,"nSyn",num2str(NoSyn),"VAF95"),[r;aa;p])
          

 
 
    
    %corr(synergyOutLeft(6).W,synergyOutRight(6).W)
end
%[r,p] = corr(wss(:,:,1), wss(:,:,2))