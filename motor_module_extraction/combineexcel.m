%source_dir = 'C:\Users\rahimig.mohammad\OneDrive - University of Florida\research\Allen works-simulations\aim1\plot synergies\only overal VAF'
%dest_dir = 'C:\Users\rahimig.mohammad\OneDrive - University of Florida\research\Allen works-simulations\aim1\plot synergies\only overal VAF'
source_files = dir(fullfile('*.xls'));
totalfile=[];
 for i = 1:length(source_files)   
     data1=string(importdata(source_files(i).name));
     for j=1:4
     ll=length(data1(:,j));
     row(1:ll,4*i-4+j) =  data1(:,j);
     end
 end
%xlswrite(strcat(ffname),row)

for k=1:54
    allVAF97(:,k)=row(:,4*k-3);
end
for k=1:54
    allVAF95(:,k)=row(:,4*k-1);
end
j=1
m=1
for i=1:54
    if contains(allVAF97(1,i),"_l")
        allVAF97_L(:,j)=allVAF97(:,i);
        j=j+1;
    else
        allVAF97_R(:,m)=allVAF97(:,i);
        m=m+1;
    end
end
j=1
m=1
for i=1:54
    if contains(allVAF95(1,i),"_l")
        allVAF95_L(:,j)=allVAF95(:,i);
        j=j+1;
    else
        allVAF95_R(:,m)=allVAF95(:,i);
        m=m+1;
    end
end
xlswrite("all_bad_muscles",row)
xlswrite("all_bad_muscles_97",allVAF97)
xlswrite("all_bad_muscles_95",allVAF95)
xlswrite("all_bad_muscles_97L",allVAF97_L)
xlswrite("all_bad_muscles_97R",allVAF97_R)
xlswrite("all_bad_muscles_95L",allVAF95_L)
xlswrite("all_bad_muscles_95R",allVAF95_R)