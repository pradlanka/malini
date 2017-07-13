Filepath=pwd;
cd(fullfile(folder_name,'output\timeseriesextract\data_save'));
Groupdir= dir;
for classes =1:length(Groupdir)-2
   %mkdir(fullfile(Filepath,'output','subjects_together',Groupdir(classes+2).name))
    d=Groupdir(classes+2,1).name;
    mkdir(fullfile(folder_name,'output','subjects_together',Groupdir(classes+2).name))
cd(Groupdir(classes+2).name)
sdir=dir;
    for k=1:(length(sdir)-2)
   %d=sdir(k+2,1).name; %directory of the current subject
    %cd(sdir(k+2,1).name);
%data_d=[];
sadi='data_';
%desicoll=strcat(sadi,d);
%strcat(sadi,d)=[];
desi=[];
TSfiles = dir('*.mat');
for l=1:length(TSfiles)
    load(TSfiles(l).name)
    desi = cat(3,desi,ts_fmri_f);
end
desicoll=strcat(sadi,d);
   cd(fullfile(folder_name,'output','subjects_together',Groupdir(classes+2).name))
    save(desicoll,'desi')
cd(fullfile(folder_name,'output','timeseriesextract','data_save',Groupdir(classes+2).name))
    %clear Dir a c d i j ts_fmri_f
    end
    cd ..
end
    
cd ..
cd ..
cd(Filepath);

