
% Static FC
%static FC is calculated using Pearson's correlation, it is a matrix of #
%of regions by # of regions by # of subjects
% n2 is # of regions, n3 is # of subjects,data is
 %extracted time series from regions, and of dimension # of time by # of regions by
 %# of subjects
 addpath(genpath(pwd));
Filepath=pwd;
cd(fullfile(folder_name,'output','subjects_together_decon'));
Groupdir= dir;
for classes =1:length(Groupdir)-2
   %mkdir(fullfile(Filepath,'output','subjects_together',Groupdir(classes+2).name))
    d=Groupdir(classes+2,1).name;
    mkdir(fullfile(folder_name,'output','deconvolution - model','sfc',Groupdir(classes+2).name))
cd(Groupdir(classes+2).name)
sdir=dir;
namem=strcat('data_decon_',d);
load(namem);
    [~,n2,n3]=size(desi_decon);
staticFC=zeros(n2,n2,n3);
for k=1:n3  
            staticFC(:,:, k) = corrcoef(  desi_decon(: ,:,k)); %using Pearson's correlation
end
%mkdir(fullfile(Filepath,'output','sfc','ptsd'))
cd(fullfile(folder_name,'output','deconvolution - model','sfc',Groupdir(classes+2).name))
bleh=strcat('staticFC_decon_',d);
save(bleh,'staticFC');
X{classes}=staticFC;
clearvars desi;
cd(fullfile(folder_name,'output','subjects_together_decon'))
end
cd(Filepath);