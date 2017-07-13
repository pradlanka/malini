addpath(genpath(pwd))
Filepath=pwd;
cd(fullfile(folder_name,'output','subjects_together'));
Groupdir= dir;
for classes =1:length(Groupdir)-2
   %mkdir(fullfile(Filepath,'output','subjects_together',Groupdir(classes+2).name))
    d=Groupdir(classes+2,1).name;
    mkdir(fullfile(folder_name,'output','sec',Groupdir(classes+2).name))
cd(Groupdir(classes+2).name)
sdir=dir;
namem=strcat('data_',d);
load(namem);
bic_smooth=find_order(desi)
order=bic_smooth;
bleh=strcat('conmat_',d);
 bleh=main_cpgc_func(desi,order)
 cd(fullfile(folder_name,'output','sec',Groupdir(classes+2).name))
sled=strcat('conmat_',d);
save(sled,'bleh');
X{classes}=bleh;
clearvars desi;
cd(fullfile(folder_name,'output','subjects_together'))
end
cd(Filepath);
 