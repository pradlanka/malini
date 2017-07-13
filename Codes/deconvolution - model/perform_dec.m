addpath(genpath(pwd));
 format compact;
Filepath=pwd;
%a=load_nii('BrainMask_05_61x73x61.img');
cc200=load_nii(templa);

cc200.img = (cc200.img>0);
b=cc200.img;
[n1,n2,n3]=size(b);
clear cc200
for j=1:n2
    for k=1:n1
        s1(j,k)=sum(sum(b(:,j,k))); %used to calculate percentage complete
    end
end
ss1=sum(sum(s1(1:n2,1:n1))); %used to calculate percentage complete

open_list=cell(46,1);  
cc200=load_nii(templa);
temp = unique(cc200.img);
bj=cell(size(temp,1)-2,1);
for j=1:(size(temp)-2)
    bj{j}=zeros(size(cc200.img));
    index = find(cc200.img==temp(j+1)) ; 
    bj{j}(index)=1;
end

cd(fullfile(folder_name,'output','deconvolution - model','run_deconv'))
Groupdir= dir;
for classes =1:length(Groupdir)-2
cd(Groupdir(classes+2).name)
sdir=dir;
for i=2:(length(sdir)-1)
    d=sdir(i+1,1).name;
    load([sdir(i+1,1).name]); %loads variable 'c' from 4-D MAT file
    timepoint=size(c,1);
    c1=zeros(timepoint,n1,n2,n3);
    H1=zeros(50,n1,n2,n3);
    ag1=zeros(n1,n2,n3);
    P1=zeros(3,n1,n2,n3);
    cd ..

    clk=clock;
    cl1=clk(4:6);
    disp(sprintf('deconvolution start time = %uh : %um : %us',cl1(1),cl1(2),round(cl1(3))))
    for j=1:n2%109
        sm = sum(sum(s1(1:j-1,1:91)));
        for k=1:n1%91
            pcent=100*(sm+sum(s1(j,1:k)))/ss1; %calculation of percentage complete
            disp(sprintf('subject(%d), voxel (%d,%d) done, percent complete = %.2f', i-1,j,k,pcent)) %display progress
            [ed H1 ag1 P1] = exec_deconv(c(:,:,j,k));
            [c1(:,:,j,k)] = single(ed);
            [HRF(:,:,j,k)] = single(H1);
            [adjust_global(:,j,k)] = single(ag1);
            [PARA(:,:,j,k)] = single(P1);
            clear ed H1 ag1 P1
                    end
    end
    clear c j k

    mkdir(fullfile(folder_name,'output','deconvolution - model','save_deconv',Groupdir(classes+2).name))
    cd(fullfile(folder_name,'output','deconvolution - model','save_deconv',Groupdir(classes+2).name))
    save(d, 'c1', '-v7.3')
    cd ..
    cd ..
   mkdir(fullfile(folder_name,'output','deconvolution - model','params',Groupdir(classes+2).name))
    cd(fullfile(folder_name,'output','deconvolution - model','params',Groupdir(classes+2).name))
    savefile=sprintf('param_%s', d);
    save(savefile, 'HRF', 'adjust_global', 'PARA', '-v7.3')
    cd ..
    
    clear HRF adjust_global PARA Dir c1 d
    clk=clock;
    cl2=clk(4:6);
    disp(sprintf('deconvolution end time = %uh : %um : %us',cl2(1),cl2(2),round(cl2(3))))

    cd(fullfile(folder_name,'output','deconvolution - model','run_deconv'))
    disp(sprintf('subject %d done', i-1))
end
end
cd ..
cd ..
cd(Filepath);
