%Run this code when everything is ready for classification
format compact; clc
%To explain briefly, the algorithm will start with an initial number of
%clusters and reduce it recursively to a final number of clusters. In each
%iteration, the classification is performed, and classification accuracy is
%obtained as well important features
hypelsvm=handles.hypelsvm;
lowv=str2double(hypelsvm{1,1});
incred=str2double(hypelsvm{2,1});
uppe=str2double(hypelsvm{3,1});
exppara=strcmp('yes',hypelsvm{4,1});
resa=handles.resa;
nocl=handles.nocl;
nofo=handles.nofo;
%below are the set of parameters
%if exist('incl')
if ~isempty(incl)
    incl=handles.incl;
init_nmc=incl;%initial number of clusters
else
    init_nmc=40;
end
if ~isempty(fincl)
    fincl=handles.fincl;
    end_nmc=fincl;
else
    end_nmc=1; %final number of clusters=1
end
if ~isempty(deccl)
    deccl=handles.deccl;
    decrease=deccl;
else
       decrease=0.5; %a decrease of 0.5 implies that number of clfincl)usters will recursively be 40->20->10->5->2
end
if ~isempty(epat)
    epat=handles.epat;
    end_paths=epat;
else
    end_paths=1;
end
topResample=resa; %no of times resampling is done.
no_classes=nocl; %insert the number of classes
kfold=nofo;
identify_important_paths = false;
Decision_Surface_Flag = DSF; % If you want the Decision Surfaces to be saved set it to 1 or 0 else.
HyperParameter_Search =HPS;  % If you want the algorithm to optimize hyperparameters set it to 1 else 0;
max_accuracy_Decision_Surface = MADC; %only saves the max acccuracy Decision Surfaces
select_significant_paths = SSP;% if you want the algoithm to perform an initial filtering step(Only significantly different paths are selected either by t-test/ANOVA on the train data)
% Input the names of  the classes
% for i=1:no_classes
%     prompt='Enter the name of the classes \n' ;
%     classname{i} = input(prompt,'s');
% end
Filepath=pwd;
%cd(folder_name);
%Groupdir=dir;
%for i=1:length(Groupdir)-2
 %   classname{(i)}=Groupdir(i+2,1).name;
%end
cd(Filepath);
%cd ..
%cd ..
%classname{1}='Controls';
%classname{2}='PTSD';
%classname{2}='EMCI';
%classname{3}='LMCI';
%classname{4}='AD';
if(arerce==1)
CurrentFolder=pwd;
addpath(genpath(fullfile(CurrentFolder,'Classifiers','LinearSVM','Classifier_Files')));
addpath(fullfile(CurrentFolder,'Classifiers','RCE','classify'));

% Ask for values hyperparmaters to be optimized using grid search

if(exppara==1)
hyperparameter(1,1)=lowv;
dummy=0;i=1;
while(dummy<=uppe)
    dummy = lowv*incred;
i=i+1;
lowv=dummy;
    if(dummy>uppe)
        break;
    end

hyperparameter(1,i+1)=dummy;
end

Hyperparameter{1}=hyperparameter;
else
Hyperparameter{1} =(lowv : incred : uppe);
end
%Hyperparameter{1} = [0.01 0.03 0.1 0.3 1 3 10 30 100];

if isempty(Hyperparameter)
    HyperParameter_Search =false;
end
Optimizable_parameters = length(Hyperparameter); % Number of Hyperparameters that can be optimized in the algorithm
mkdir(fullfile(folder_name, filesep,'Rce_Results',discos));


% S = 'Y';
% while (S == 'Y' ||  S=='y')
% if (HyperParameter_Search ==1)
%     for i=1:Optimizable_parameters
%     fprintf('The Default Regularization Parameter Grid Search\n');
%     disp( Hyperparameter{i});
%     Parameter_Change= input('Do you want to change the default search space? Type Y or y if yes, else press N or n\n','s');
%     if (Parameter_Change =='Y' || Parameter_Change =='y')
%        prompt = ('Please enter the parameter search space\n');
%         Hyperparameter{i} = input(prompt);
%     end
%     end
% end
%[~,~,input_data]=xlsread('Table.xlsx','SFC');
%[~,~,input_data]=xlsread('train.xlsx');
%load input_data.mat     %load the data contained in inputdata
v= ver;
if any(strcmp('Parallel Computing Toolbox', {v.Name}));
      num_workers= parpool('local');
    if num_workers==0
        parpool open;
    end
end
cd(fullfile(CurrentFolder,'Classifiers','RCE','classify'))        %change the current folder to classify
mkdir(fullfile(folder_name,'Rce_Results',discos,'Accuracy_vals'));
if max_accuracy_Decision_Surface 
   mkdir(fullfile(folder_name,'Rce_Results',discos,'MaxAccuracy_vals'));
end
% % Run the first run to estimate the size of the Decision Surface and the time taken for the algorithm
% InitResample=1;
% if (Decision_Surface_Flag) &&  (HyperParameter_Search)
%     mkdir(fullfile(CurrentFolder,'Results','Decision_Surfaces'))
% [time, Decision_Surf_size] = Run_SVM_RCE(CurrentFolder,input_data,init_nmc,end_nmc,decrease,end_paths,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,max_accuracy_Decision_Surface,HyperParameter_Search,InitResample,1,Hyperparameter);
% elseif(Decision_Surface_Flag) &&  ~(HyperParameter_Search)
%     mkdir(fullfile(CurrentFolder,'Results','Decision_Surfaces'))
%     [time, Decision_Surf_size] = Run_SVM_RCE(CurrentFolder,input_data,init_nmc,end_nmc,decrease,end_paths,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,max_accuracy_Decision_Surface,HyperParameter_Search,InitResample,1);
% elseif ~(Decision_Surface_Flag) &&  (HyperParameter_Search)
%     [time] = Run_SVM_RCE(CurrentFolder,input_data,init_nmc,end_nmc,decrease,end_paths,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,max_accuracy_Decision_Surface,HyperParameter_Search,InitResample,1,Hyperparameter);  
% else
%     [time] = Run_SVM_RCE(CurrentFolder,input_data,init_nmc,end_nmc,decrease,end_paths,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,max_accuracy_Decision_Surface,HyperParameter_Search,1,InitResample);
% end

% % Prompt the HardDrive size, RAM Required and the total time taken
% cd(CurrentFolder);
% fileID = fopen('Initial_Time_&_Space_estimate.txt','w');
% fprintf('The total running time for the algorithm = %5.3f hours\n',time*topResample/(num_workers*3600));
% fprintf(fileID,'The total running time for the algorithm = %5.3f hours\n',time*topResample/(num_workers*3600));
% if (Decision_Surface_Flag)
% fprintf('The total space on Harddrive required = %4.2f GB\n',Decision_Surf_size*topResample/2^30); 
% fprintf(fileID,'The total space on Harddrive required = %4.2f GB\n',Decision_Surf_size*topResample/2^30); 
% fprintf('The total RAM required = %3.2f GB\n',Decision_Surf_size/2^30);
% end 
% fclose(fileID);
% S = lower(input('Do you want to continue? if  Yes Press Y. If you want to change the settings press N or press ctrl+C to cancel the program\n','s'));
% if S=='y'
%     break;
% end
% 
% if S== 'n'
% DS = lower(input('Do you want to change the Decision_Surface_Flag?\n press Y for Yes or N for No\n','s'));
% if DS == 'y'
%     Decision_Surface_Flag = ~Decision_Surface_Flag;    
% end
% DS = lower(input('Do you want to change the HyperParameter_Search_Flag?\n press Y for Yes or N for No\n','s'));
% if DS == 'y'
%     HyperParameter_Search_Flag = ~HyperParameter_Search_Flag;    
% end
% end
% end
if (Decision_Surface_Flag)
  mkdir(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces'))
end

cd(fullfile(CurrentFolder,'Classifiers','RCE','classify'))
%Run the algorithm
if (HyperParameter_Search)
   [time] = Run_SVM_RCE(CurrentFolder,input_data,init_nmc,end_nmc,decrease,end_paths,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,max_accuracy_Decision_Surface,HyperParameter_Search,topResample,0,Hyperparameter);  
else
    [time] = Run_SVM_RCE(CurrentFolder,input_data,init_nmc,end_nmc,decrease,end_paths,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,max_accuracy_Decision_Surface,HyperParameter_Search,topResample,0);
end
cd(fullfile(folder_name,'Rce_Results',discos));
fileID = fopen('Total_time.txt','w');
fprintf(fileID,'Total time = %5.3f hours\n',time/3600);
fclose(fileID);
fprintf('Total time = %5.3f hours\n',time/3600);
cd(CurrentFolder);
%rmdir(fullfile(CurrentFolder,'Rce_Results',discos,'MaxAccuracy_vals'));
rmpath(genpath(fullfile(CurrentFolder,'Classifiers','LinearSVM','Classifier_Files')));
rmpath(fullfile(CurrentFolder,'Classifiers','RCE','classify'));
if any(strcmp('Parallel Computing Toolbox', {v.Name}));
delete(gcp('nocreate'));
end
datestr(now)
cd(CurrentFolder);

else
    Currentfolder = pwd;
addpath(genpath(fullfile(Currentfolder,'Classifiers','LinearSVM','Classifier_Files')));
addpath(fullfile(Currentfolder,'Classifiers','Parameter_Optimization','classify'));
%Hyperparameter{1} =(lowv : incred : uppe);

if(exppara==1)
hyperparameter(1,1)=lowv;
dummy=0;i=1;
while(dummy<=uppe)
    dummy = lowv*incred;
i=i+1;
lowv=dummy;
    if(dummy>uppe)
        break;
    end

hyperparameter(1,i+1)=dummy;
end

Hyperparameter{1}=hyperparameter;
else
Hyperparameter{1} =(lowv : incred : uppe);
end
Optimizable_parameters = length(Hyperparameter);
if isempty(Hyperparameter)
    HyperParameter_Search =false;
end
mkdir(fullfile(folder_name, filesep,'Outside_Rce_Results',discos));
mkdir(fullfile(folder_name, filesep,'Outside_Rce_Results',discos,filesep,'Accuracy'));

if identify_important_paths
    mkdir(fullfile(folder_name,'Outside_Rce_Results',discos,'Feature_Significance'))
end

v= ver;
if any(strcmp('Parallel Computing Toolbox', {v.Name}));
      num_workers= parpool('local');
    if num_workers==0
        parpool open;
    end
end


if (Decision_Surface_Flag)
    mkdir(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Decision_Surfaces'));
end
if (HyperParameter_Search)
    mkdir(fullfile(folder_name, filesep,'Outside_Rce_Results',discos,filesep,'bestparameters'));
end

cd(fullfile(Currentfolder,'Classifiers','Parameter_Optimization','classify'))
%Run the algorithm
if (HyperParameter_Search)
   [time] = Run_Classifier(Currentfolder,input_data,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,HyperParameter_Search,identify_important_paths,topResample,0,Hyperparameter);  
else
   [time] = Run_Classifier(Currentfolder,input_data,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,HyperParameter_Search,identify_important_paths,topResample,0);
end
cd(fullfile(Currentfolder,'Outside_Rce_Results',discos));
fileID = fopen('Total_time.txt','w');
fprintf(fileID,'Total time = %5.3f hours\n',time/3600);
fclose(fileID);
fprintf('Total time = %5.3f hours\n',time/3600);
rmpath(genpath(fullfile(Currentfolder,'Classifiers','LinearSVM','Classifier_Files')));
rmpath(fullfile(Currentfolder,'Classifiers','Parameter_Optimization','classify'));
if any(strcmp('Parallel Computing Toolbox', {v.Name}));
    delete(gcp('nocreate'));
end
datestr(now);
cd(Currentfolder);
end