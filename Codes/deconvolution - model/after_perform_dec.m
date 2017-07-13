addpath(genpath(pwd));
Filepath = pwd;
cd(fullfile(folder_name,'output\deconvolution - model\save_deconv'));
Groupdir= dir;
for classes =1:length(Groupdir)-2
cd(Groupdir(classes+2).name)
sdir=dir;
   for k=1:(length(sdir)-2)
    d=sdir(k+2).name; %directory of the current subject
   % cd(sdir(k+2,1).name); %goes into the directory of the current subject
   % Create a mat filename, and load it into a structure called matData.
	matFileName = sprintf('Sub_00%d.mat', k);
	if exist(matFileName, 'file')
		load(matFileName);
    c1= permute(c1,[2,3,4,1]);
    cd(data_folder_name);
    Groupdir=dir;
    cd(Groupdir(3,1).name);
    nadir=dir;
    cd(nadir(3,1).name);
    fakir=dir;
    nig=fakir(3,1).name;
    sis=load_nii(nig);
    %sis=load_nii('Filtered_4DVolume.nii');
   sis.img=c1;
    end
    cd ..
    cd ..
    mkdir(fullfile(folder_name,'output','deconvolution - model','post_decon',Groupdir(classes+2).name))
   cd(fullfile(folder_name,'output','deconvolution - model','post_decon',Groupdir(classes+2).name))
   save_nii(sis,strcat('postdecon','_',sdir(k+2).name,'.nii'));
   cd ..
   cd ..
   cd(fullfile(folder_name,'output','deconvolution - model','save_deconv',Groupdir(classes+2).name))
  clear Dir a c d i j 
    disp(sprintf('Group %d:subject %d done', classes,k))
   end
       cd ..
end
cd ..
cd ..
cd(Filepath);
