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
CI=95
for n = 1:length(source_files);
    filename=source_files(n);
    datan(1,n)=filename;
    solution=importdata(strcat(filename,".sto"));
    syndata=importdata(strcat(filename,"_19-Jul-2023.mat"));
    leftsyn=syndata.synergyOutLeft;
    Rightsyn=syndata.synergyOutRight;
DATA = Rightsyn;
% how many replicates?
if ~isfield(DATA,'sortedVAF_original')
    nreps = 100;
else
    nreps = size(DATA(1).sortedVAFmus_original,2);
end
temp = (100 - CI)/(2*100);
up = ceil((1 - temp)*nreps);
bot = floor(temp*nreps);

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



errorbar(1:length(DATA),[DATA(:).VAF], ([DATA(:).VAF]-lb), (ub-[DATA(:).VAF]));
hold on
grid on
%errorbar(1:length(DATA), x2, (x2-lb2), (ub2-x2),'color','r');
xx=[0 10];
yy= [95 95];
%line(xx,yy,'Color','red','LineStyle','--');
xlim([3 7]);

title([strcat("@",' Overal VAF with') num2str(CI) '% CI']);

end
legend(source_files);