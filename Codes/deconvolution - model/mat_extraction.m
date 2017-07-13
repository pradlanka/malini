% Extract 4-D mat file (timex91x109x91) from 4-D NIFTI images
addpath(genpath(pwd));
Filepath=pwd;
format compact;
cc200=load_nii(templa);
cc200.img = (cc200.img>0);
temp=(cc200.img~=0);
cd(data_folder_name);
Groupdir= dir;
for classes =1:length(Groupdir)-2
cd(Groupdir(classes+2).name)
sdir=dir;
for i=1:(length(sdir)-2)
    d=sdir(i+2,1).name;
    cd(sdir(i+2,1).name);
    Dir=dir;
    a=load_nii([Dir(3).name]);
    [x,y,z,timepoints]=size(a.img);
    c=zeros(timepoints,x ,y,z);
    for j=1:timepoints
        c(j,:,:,:)=a.img(:,:,:,j).*single(cc200.img);
        sprintf('subject %d, timepoint %d',i,j)
    end
    clear a
    cd ..
    cd ..
    cd ..
    
    
    mkdir(fullfile(folder_name,'output','deconvolution - model','run_deconv',Groupdir(classes+2).name))
    cd(fullfile(folder_name,'output','deconvolution - model','run_deconv',Groupdir(classes+2).name))
    save(d, 'c', '-v7.3')
    cd ..
    cd ..    
    cd(fullfile(fullfile(data_folder_name,Groupdir(classes+2).name)))
    clear Dir a c d j
    disp(sprintf('Group %d:subject %d done', classes,i))
end
cd ..

end
cd ..
cd ..
cd(Filepath);
%%