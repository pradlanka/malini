% Extract 190 channel time series from Tx91x109x91 fMRI images (T=number of timepoints)
 format compact;
addpath(genpath(pwd));
Filepath = pwd;
cc200=load_nii(templa);
temp=unique(cc200.img);
b=cell(size(temp,1)-1,1);
for j=1:(size(temp)-1)
    b{j}=zeros(91,109,91);
    index = find(cc200.img==temp(j+1));
    b{j}(index)=1;
    clear index
end
cd(fullfile(folder_name,'output\deconvolution - model\post_decon'));
Groupdir= dir;
for classes =1:length(Groupdir)-2
 mkdir(fullfile(folder_name,'output','deconvolution - model','post_decon_extract',Groupdir(classes+2).name))
    cd(Groupdir(classes+2).name)
sdir=dir;
for k=1:(length(sdir)-2)
    %d=sdir(k+2,1).name; %directory of the current subject
    %cd(sdir(k+2,1).name); %goes into the directory of the current subject
    fMRI=load_nii(sdir(k+2).name);
    timepoint=size(fMRI.img,4); %For ex., if there are 1000 img and 1000 hdr files then there will be a total of 2002 files, of which 2 of them refer to prev directory and root directory
    p_ts_fmri_f=zeros(timepoint, size(temp,1)-1); %averaged timeseries
    c=zeros(timepoint, size(cc200.img,1),size(cc200.img,2),size(cc200.img,3) ); %#ok<NASGU> % initialized: Tx91x109x91 matrix
      Volfile = fMRI.img;
    for i=1:timepoint
        for j=1:(size(temp)-1) % j=1:190 here
            c=Volfile(:,:,:,i).*single(b{j}); %selecting only those voxels defined by the cc200 region b{j} (these will be 1 and other voxels will be 0)
            index= find(b{j}~=0); %3D coordinates of the selected region
            p_ts_fmri_f(i,j) = mean(c(index)); %averaged timeseries = mean of all voxel values in that region
            clear c
        end % the process is repeated for all 190 channels
        disp(sprintf('subject %d, timepoint %d',k,i))
    end % the process is repeated for all timepoints
    cd ..
    cd ..
   Filename=strsplit(sdir(k+2).name,'.');
   cd(fullfile(folder_name,'output','deconvolution - model','post_decon_extract',Groupdir(classes+2).name))
    save(strcat(Filename{1},'.mat'),'p_ts_fmri_f');
    %save(d, 'p_ts_fmri_f')
    
    cd ..
    cd(fullfile(folder_name,'output','deconvolution - model','post_decon',Groupdir(classes+2).name))
    clear Dir a c d i j p_ts_fmri_f
    disp(sprintf('group %d:subject %d done', classes,k))
end
cd ..

end
cd ..
cd ..
cd(Filepath);


%%