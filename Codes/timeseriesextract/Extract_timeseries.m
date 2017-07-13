% Extract 190 channel time series from Tx91x109x91 fMRI images (T=number of timepoints)
format compact;
addpath(genpath(pwd));
Filepath = pwd;
%cc200=load_nii('cc200_91x109x91.nii');
cc200=load_nii(templa);
%cc200=load_nii('KKI_0050822_rois_ho.1D');
szei=size(cc200.img);
temp=unique(cc200.img);
b=cell(size(temp,1)-1,1);

for j=1:(size(temp)-1)
    % b{j}=zeros(91,109,91);
    b{j}=zeros(szei);
    index = find(cc200.img==temp(j+1));
    b{j}(index)=1;
    clear index
end
%cd(fullfile(folder_name,'Preprocessed-data'));
cd(data_folder_name);
Groupdir= dir;
for classes =1:length(Groupdir)-2
     mkdir(fullfile(folder_name,'output','timeseriesextract','data_save',Groupdir(classes+2).name));
cd(Groupdir(classes+2).name)
sdir=dir;
for k=1:(length(sdir)-2)
    d=sdir(k+2,1).name; %directory of the current subject
    cd(sdir(k+2,1).name); %goes into the directory of the current subject
    niidir=dir;
     fMRI=load_nii(niidir(3,1).name);
    %fMRI=load_nii('Filtered_4DVolume.nii');
    timepoint=size(fMRI.img,4); %For ex., if there are 1000 img and 1000 hdr files then there will be a total of 2002 files, of which 2 of them refer to prev directory and root directory
    ts_fmri_f=zeros(timepoint, size(temp,1)-1); %averaged timeseries
    c=zeros(timepoint, size(cc200.img,1),size(cc200.img,2),size(cc200.img,3) ); %#ok<NASGU> % initialized: Tx91x109x91 matrix
      Volfile = fMRI.img;
    for i=1:timepoint
        for j=1:(size(temp)-1) % j=1:190 here
            c=Volfile(:,:,:,i).*single(b{j}); %selecting only those voxels defined by the cc200 region b{j} (these will be 1 and other voxels will be 0)
            index= find(b{j}~=0); %3D coordinates of the selected region
            ts_fmri_f(i,j) = mean(c(index)); %averaged timeseries = mean of all voxel values in that region
            clear c
        end % the process is repeated for all 190 channels
        disp(sprintf('subject %d, timepoint %d',k,i))
    end % the process is repeated for all timepoints
    cd ..
    cd ..
    cd ..
    %cd(fullfile(Filepath,'Codes','timeseriesextract','data_save',Groupdir(classes+2).name))
    cd(fullfile(folder_name,'output','timeseriesextract','data_save',Groupdir(classes+2).name))
    save(d, 'ts_fmri_f')
    cd ..
    cd ..
    cd(fullfile(data_folder_name,Groupdir(classes+2).name));
    clear Dir a c d i j ts_fmri_f
    disp(sprintf('Group %d:subject %d done', classes,k));
    
end
cd ..
%movefile(ef,'Preprocessed-data');
cd(data_folder_name);
end
cd ..
cd ..
cd(Filepath);


%%