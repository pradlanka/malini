addpath(genpath(pwd))
Filepath=pwd;
cd(fullfile(folder_name,'output','subjects_together_decon'));
Groupdir= dir;
for classes =1:length(Groupdir)-2
   %mkdir(fullfile(Filepath,'output','subjects_together',Groupdir(classes+2).name))
    d=Groupdir(classes+2,1).name;
    mkdir(fullfile(folder_name,'output','deconvolution - model','dfc',Groupdir(classes+2).name))
cd(Groupdir(classes+2).name)
sdir=dir;
namem=strcat('data_decon_',d);
load(namem);
bleh=strcat('dynamicFC_decon_',d);
 n4=10;n5=40;
 %[smooth_window_length_Controls, dynamicFC_Controls]=adf_slidingwindowFC(data__dfc_Controls,n4,n5)
 [smooth_window_length_Controls, bleh]=adf_slidingwindowFC(desi_decon,n4,n5)
 bleh=permute(bleh,[2,3,4,1]);
 bleh=var(bleh,[],4);
 cd(fullfile(folder_name,'output','deconvolution - model','dfc',Groupdir(classes+2).name))
sled=strcat('dynamicFC_decon_',d);
save(sled,'bleh');
X{classes}=bleh;
clearvars desi_decon;
cd(fullfile(folder_name,'output','subjects_together_decon'))
end
cd(Filepath);
