%% plotting w and reconstructed Cs based on the speed and R of Cs when compared to the reference
allCrlsqnnM=importdata("all_recons_C_r_finalJan.xls");
allCllsqnnM=importdata("all_recons_C_l_finalJan.xls")
syndata1=importdata('11000normalsecond_19-Jul-2023.mat');
NoSyn=4
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
muscleNames1={'piri_r','psoas_r','sart_r','tfl_r','iliacus_r','addbrev_r','addlong_r','addmagIsch_r','addmagDist_r','addmagMid_r','addmagProx_r',...
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

wr=syndata1.synergyOutRight(4).W;
wl=syndata1.synergyOutLeft(4).W;

for j=1:2
    if j==1
        ws=wr;
        cs=allCrlsqnnM;
        xx=syndata1.synergyOutRight(4).C;
        side="_R"
    elseif j==2
        ws=wl;
        cs=allCllsqnnM;
        xx=syndata1.synergyOutLeft(4).C;
        side="_L"
    end
            figg(j)=figure('Renderer', 'painters', 'Position', [10 10 2000 1200]);
            newfname=strcat("WCR")
            for i= 1:NoSyn
                subplot(NoSyn,6,[6*i-5 6*i-4])
                %hold on; box on;
                X=[1:43];
                Y=ws(:,i);
                bar(X,Y,'FaceColor','#3771C8')
                ylabel('Activation','FontWeight','b')
                set(gca,'xtick',1:length(muscleNames1))
                set(gca,'xticklabel',muscleNames1(:),'FontSize',7)
                xtickangle(90)
                title(strcat("Ws composition",newfname,side));
                
                %
                subplot(NoSyn,6,[6*i-3 6*i-1])
                linest=["-","--",":"];
                lineWd=[0.8,0.8,1.2]
                for jjj=1:3
                %plot(time(:)/time(end)*100,cs(4*9*(jjj-1)+i,:),linest(jjj),'LineWidth',lineWd(jjj)) %% for across speeds normal
                %plot(time(:)/time(end)*100,cs(4*(jjj-1)+i+18*4,:),linest(jjj),'LineWidth',lineWd(jjj)) %% for across step width normal
                plot(time(:)/time(end)*100,cs(18*4+4*3*(jjj-1)+i,:),linest(jjj),'LineWidth',lineWd(jjj)) %% for across asymmetry
                hold on
                end
                hold off
                xlabel('Stride (%)','FontWeight','b')
                ylabel('Activation','FontWeight','b')
                % colorful colors
                %newcolors = {'#DC143C','#9932CC','#FF8C00','#FFD700','#A0522D','#008000','#32CD32','#008B8B','#1E90FF','#000080'};
                %shades of one color
                newcolors = {'#3771C8','#3771C8','#3771C8'};
                colororder(newcolors)
                legend("0%","15%","30%")
                legend('Location','eastoutside')
                
                
                %ylim([-30 30])
                title(strcat("Cs activations for 1.45m/s",newfname,side));
                subplot(NoSyn,6,6*i)
                %hold on; box on;
                X=[1:3];
                Y=[];
                for jjj=1:3
                %[CREL,CPEL] = corr(xx(i,:)',cs(4*9*(jjj-1)+i,:)');
                %[CREL,CPEL] = corr(xx(i,:)',cs(4*(jjj-1)+i+18*4,:)');
                [CREL,CPEL] = corr(xx(i,:)',cs(18*4+4*3*(jjj-1)+i,:)');
                Y=[Y,CREL];
                end
                labels=["0.1m","0.2m","0.3m"];
                h=bar(X,Y,'FaceColor','	#0000CD')
                ylabel('R value','FontWeight','b')
                set(gca,'xtick',1:3)
                set(gca,'xticklabel',labels)
                minBar = min(Y);
                maxBar = max(Y);
                ylim([0.95*minBar, 1.05*maxBar])
                xtickangle(90)
                title(strcat("R of Cs",newfname,side))
                colororder(newcolors)
            end
                %exportgraphics(figg(j),strcat("WandCs","speedcomp","nSyn",num2str(NoSyn),side,'.jpg'),'Resolution',500)
                
            end
        
    


