format compact; clc;
%%% Initial Parameters
resa=handles.resa;
nocl=handles.nocl;
nofo=handles.nofo;
hyperotfo=handles.hyperotfo;
lowv=str2double(hyperotfo{1,1});
incred=str2double(hyperotfo{2,1});
uppe=str2double(hyperotfo{3,1});
exppara=strcmp('yes',hyperotfo{4,1});
topResample=resa; %no of times resampling is done.
no_classes=nocl; %insert the number of classes
kfold=nofo;
Decision_Surface_Flag = DSF; % If you want the Decision Surfaces to be saved set it to 1 or 0 else.
HyperParameter_Search =HPS;  % If you want the algorithm to optimize hyperparameters set it to 1 else 0;
select_significant_paths = SSP; % If you want the algoithm to perform an initial filtering step(Only significantly different paths are selected)
identify_important_paths = false;% Input the names of  the classes
% for i=1:no_classes
%     prompt='Enter the name of the classes \n' ;
%     classname{i} = input(prompt,'s');
% end

if(arerce==0)

%%% Add path for the toolbox classifier
CurrentFolder = pwd;
addpath(genpath(fullfile(CurrentFolder,'Classifiers','Rotation Forest','Classifier_Files')));
addpath(fullfile(CurrentFolder,'Classifiers','Parameter_Optimization','classify'));
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
%Hyperparameter{1} = 1:4:40;

Optimizable_parameters = length(Hyperparameter);
if isempty(Hyperparameter)
    HyperParameter_Search =false;
end
mkdir(fullfile(folder_name, filesep,'Outside_Rce_Results',discos));
mkdir(fullfile(folder_name, filesep,'Outside_Rce_Results',discos,filesep,'Accuracy'));

if identify_important_paths
    mkdir(fullfile(folder_name,'Outside_Rce_Results',discos,'Feature_Significance'))
end

%load input_data.mat     %load the data contained in inputdata
v= ver;
if any(strcmp('Parallel Computing Toolbox', {v.Name}));
      num_workers= parpool('local');
    if num_workers==0
        parpool open;
    end
end
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


% cd(fullfile(CurrentFolder,'Parameter_Optimization','classify'))        %change the current folder to classify
% % Run the first run to estimate the size of the Decision Surface and the time taken for the algorithm
% InitResample=1;
% if (Decision_Surface_Flag && HyperParameter_Search)
%     mkdir(fullfile(CurrentFolder,filesep,'Results',filesep,'Decision_Surfaces'));
%     mkdir(fullfile(CurrentFolder, filesep,'Results',filesep,'bestparameters'));
% [time, Decision_Surf_size] = Run_Classifier(CurrentFolder,input_data,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,HyperParameter_Search,identify_important_paths,InitResample,1,Hyperparameter);
% elseif(Decision_Surface_Flag && ~HyperParameter_Search)
%     mkdir(fullfile(CurrentFolder,filesep,'Results',filesep,'Decision_Surfaces'));
%     [time, Decision_Surf_size] = Run_Classifier(CurrentFolder,input_data,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,HyperParameter_Search,identify_important_paths,InitResample,1);
% elseif(~Decision_Surface_Flag && HyperParameter_Search)
%     mkdir(fullfile(CurrentFolder, filesep,'Results',filesep,'bestparameters'));
%     [time] = Run_Classifier(CurrentFolder,input_data,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,HyperParameter_Search,identify_important_paths,InitResample,1,Hyperparameter);  
% else
%     [time] = Run_Classifier(CurrentFolder,input_data,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,HyperParameter_Search,identify_important_paths,1,InitResample);
% end
% 
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
% S = input('Do you want to continue? if  Yes Press Y. If you want to change the settings press N or press ctrl+C to cancel the program\n','s');
% if S== 'Y' || S== 'y'
%     break;
% end
% 
% if S== 'N' || S== 'n'
% DS = input('Do you want to change the Decision_Surface_Flag?\n press Y for Yes or N for No\n','s');
% if DS == 'Y' ||DS == 'y'
%     Decision_Surface_Flag = ~Decision_Surface_Flag;    
% end
% DS = input('Do you want to change the HyperParameter_Search_Flag?\n press Y for Yes or N for No\n','s');
% if DS == 'Y' || DS =='y'
%     HyperParameter_Search_Flag = ~HyperParameter_Search_Flag;    
% end
% end
% end

if (Decision_Surface_Flag)
    mkdir(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Decision_Surfaces'));
end
if (HyperParameter_Search)
    mkdir(fullfile(folder_name, filesep,'Outside_Rce_Results',discos,filesep,'bestparameters'));
end

cd(fullfile(CurrentFolder,'Classifiers','Parameter_Optimization','classify'))
%Run the algorithm
if (HyperParameter_Search)
   [time] = Run_Classifier(CurrentFolder,input_data,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,HyperParameter_Search,identify_important_paths,topResample,0,Hyperparameter);  
else
   [time] = Run_Classifier(CurrentFolder,input_data,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,HyperParameter_Search,identify_important_paths,topResample,0);
end
cd(fullfile(folder_name,'Outside_Rce_Results',discos));
fileID = fopen('Total_time.txt','w');
fprintf(fileID,'Total time = %5.3f hours\n',time/3600);
fclose(fileID);
fprintf('Total time = %5.3f hours\n',time/3600);
rmpath(genpath(fullfile(CurrentFolder,'Classifiers','Rotation Forest','Classifier_Files')));
rmpath(fullfile(CurrentFolder,'Classifiers','Parameter_Optimization','classify'));
if any(strcmp('Parallel Computing Toolbox', {v.Name}));
    delete(gcp('nocreate'));
end
datestr(now);
cd(CurrentFolder)

else
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
    CurrentFolder=pwd;
addpath(genpath(fullfile(CurrentFolder,'Classifiers','Rotation Forest','Classifier_Files')));
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
if isempty(Hyperparameter)
    HyperParameter_Search =false;
end
Optimizable_parameters = length(Hyperparameter); % Number of Hyperparameters that can be optimized in the algorithm
mkdir(fullfile(folder_name, filesep,'Rce_Results',discos));

v= ver;
if any(strcmp('Parallel Computing Toolbox', {v.Name}));
      num_workers= parpool('local');
    if num_workers==0
        parpool open;
    end
end
cd(fullfile(CurrentFolder,'Classifiers','RCE','classify'))        %change the current folder to classify
mkdir(fullfile(folder_name,'Rce_Results',discos,'Accuracy_vals'))
if max_accuracy_Decision_Surface 
   mkdir(fullfile(folder_name,'Rce_Results',discos,'MaxAccuracy_vals'))
end

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
rmpath(genpath(fullfile(CurrentFolder,'Classifiers','Rotation Forest','Classifier_Files')));
rmpath(fullfile(CurrentFolder,'Classifiers','RCE','classify'));
if any(strcmp('Parallel Computing Toolbox', {v.Name}));
    delete(gcp('nocreate'));
end
datestr(now)
cd(CurrentFolder);
end