Filepath=pwd;
cd(fullfile(folder_name,'output\deconvolution - model\post_decon_extract'));
Groupdir= dir;
for classes =1:length(Groupdir)-2
   %mkdir(fullfile(Filepath,'output','subjects_together',Groupdir(classes+2).name))
    d=Groupdir(classes+2,1).name;
    mkdir(fullfile(folder_name,'output','subjects_together_decon',Groupdir(classes+2).name))
cd(Groupdir(classes+2).name)
sdir=dir;
    for k=1:(length(sdir)-2)
   %d=sdir(k+2,1).name; %directory of the current subject
    %cd(sdir(k+2,1).name);
%data_d=[];
sadi='data_decon_';
%desicoll=strcat(sadi,d);
%strcat(sadi,d)=[];
desi_decon=[];
TSfiles = dir('*.mat');
for l=1:length(TSfiles)
    load(TSfiles(l).name)
   desi_decon = cat(3,desi_decon,p_ts_fmri_f);
end
desicoll_decon=strcat(sadi,d);
   cd(fullfile(folder_name,'output','subjects_together_decon',Groupdir(classes+2).name))
    save(desicoll_decon,'desi_decon')
cd(fullfile(folder_name,'output','deconvolution - model','post_decon_extract',Groupdir(classes+2).name))
    %clear Dir a c d i j ts_fmri_f
    end
    cd ..
end
    
cd ..
cd ..
cd(Filepath);
