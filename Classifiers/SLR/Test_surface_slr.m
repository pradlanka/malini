if(arerce==0)
Currentfolder=pwd;
addpath(genpath(fullfile(Currentfolder,'Classifiers','SLR','Classifier_Files')));

% Load a few variables from  Results
global discos folder_name;
cd(fullfile(folder_name,'Outside_Rce_Results',discos));
load Results.mat
no_iter =kfold*topResample;
%clearvars -except classes classname no_iter paths accuracy_final FinalConfmat topResample select_significant_paths
%clearvars -except classes classname no_iter paths FinalConfmat topResample select_significant_paths
cd(Currentfolder);
%Currentfolder=pwd;

% Load test data
%cd ..
%load test_data.mat 
%test_dat=handles.test_dat;
%test_data=test_dat;
%[~,~,test_data]=xlsread('train.xlsx');
cd(Currentfolder);
test_data=handles.test_data;
test_data_targets = (test_data(2,3:end)); 
coordinates = (test_data(3:end,2));
test_data= test_data(3:end,3:end);
test_data_inp =(cell2mat(test_data))';

% Convert group names to group nos for the testing data 
no_classes=size(classname,2);  
group_no_targets=zeros(size(test_data_targets)) ;
for ns=1:1:length(test_data_targets);
    for ng=1:no_classes
        if(strcmp(test_data_targets(ns),classname{ng})==1)
        group_no_targets(ns)=ng;
        end
    end
end
 cp = classperf(group_no_targets);
% Predict the values of  Test data
classes_predicted=[]; test_accuracy_ite =[];
cd(fullfile(folder_name,'Outside_Rce_Results',discos,'Decision_Surfaces'))
h = waitbar(0,'Initializing waitbar...');
list = dir('Decision_Surface*.mat');
  for k=1:length(list)
        load(list(k).name)  
        class_act = zeros(size(Decision_Surface,1),size(test_data_inp,1));
        test_accuracy_resample = zeros(size(Decision_Surface,1),1);
     for i=1:size(Decision_Surface,1)
               model = Decision_Surface{i};
               if select_significant_paths
               index = ismember(coordinates,Paths{i});
               test_data_red =test_data_inp(:,index);
               else
               test_data_red =test_data_inp;
               end
              [class_act(i,:)] = classifier_predict(model,test_data_red);
                 test_accuracy_resample(i) =mean(class_act(i,:)==group_no_targets);
      end
     waitbar(k/topResample,h,sprintf('Calculating... %d%%',k*100/topResample));
          clear Decision_Surface
         classes_predicted = cat(1,classes_predicted, class_act);
         test_accuracy_ite = cat(1,test_accuracy_ite,test_accuracy_resample); 
  end

test_data_predict = mode(classes_predicted,1);
classperf(cp,test_data_predict);
test_data_accuracy_actual = mean(test_accuracy_ite);
test_data_accuracy_std = std(test_accuracy_ite);
test_data_accuracy_voting = cp.CorrectRate;
test_data_specificity_voting = cp.Specificity;
test_data_sensitivity_voting =cp.Sensitivity;
Confmat_test=confusionmat(group_no_targets,test_data_predict,'order',classes);
clear Decision*
[train_data_acc] = mean(accuracy_final,1);
train_data_std = std(accuracy_final);
Conf_mat_train = FinalConfmat;

for j=1:no_classes
    temp(2*(j-1)+1) = Conf_mat_train(j,j)/sum(squeeze(Conf_mat_train(j,:)));
            temp(2*j) = Confmat_test(j,j)/sum(squeeze(Confmat_test(j,:)));
    assignin('base',strcat('acc_train_',classname{j}),Conf_mat_train(j,j)/sum(squeeze(Conf_mat_train(j,:))));
    assignin('base',strcat('acc_test_',classname{j}), Confmat_test(j,j)/sum(squeeze(Confmat_test(j,:))));
end
balanced_CV_accuracy=mean(temp(1:2:end));
balanced_test_accuracy=mean(temp(2:2:end));
rmpath(genpath(fullfile(Currentfolder,'Classifiers','SLR','Classifier_Files')));
cd(fullfile(folder_name,'Outside_Rce_Results',discos));
close(h);
save Test_results -regexp ^(?!(CurrentFolder|ltinput_data|grpdat|edc|hope|hopet|input_data|test_data|cdec|grpdats|grpdatind|X_train|X_test|test_dat|kfold|topResample|classes|classname|no_iter|RCE_steps|parameter_optimized_Decision_Surface|max_accuracy_Decision_Surface|use_accuracy|acc_for_each_iteration|acc_for_each_iteration_opt|FinalConfmat|FinalConfmat_opt|time|timer|tmse|TOP_RESAMPLE|total_freq|ts_fmri_f|TSfiles|uppe|v|X1|X2|y1|y2|YourFile|clas2|clas3|clas4|clas5|clas6|clas7|clas8|clas9|clas10|clas11|clas12|clas13|clas14|clas15|clas16|clas17|clas18|clas19|clas20|clas21|classes_predicted|data|Data|data_test|data_train|data_train_tmp|Decision_name|Decision_surf_ite|Decision_surf_opt_ite|decrease|desi|desicoll|dfg|DSF|end_nmc|end_paths|Estimate_size_time|eventdata|fid|fileID|Filepath|first|folder_name|fout|fpaths|freq|fsamples|gfd|gindex|group_no_targets|Groupdir|groups_ID|groups_name|groups_test|groups_train|grp|gsorted|h|handles|hObject|hps|HPS|Hyperparameter|Hyperparameter_Search|i|idod|iip|IIP|incred|indices2|indx|init_nmc|input_data|j|jj|k|kk|l|Labels|last_num_clusters|list|loc|lowv|madc|MADC|lowto|incredto|uppeto|n|n2|n3|n_samples|ng|nmc|no_classes|no_clusters|nocl|nofo|ns|resa|resample|row|sadi|score|sdir|select_significant_paths|selected_paths_resample|significant_paths|specificity|specificity_opt|splirat|split_point|Spot_id|Spot_ID|Spot_ID_acc|Spot_ID_acc_opt|Spot_ID_ite|Spot_ID_perf|Spot_ID_TF|Spot_ID_tmp|Spotid|ssp|SSP|st|staticFC_Controls|staticFC_PTSD|stfc|sth|Table|temp|test)$).
cd(Currentfolder);
else
   
    %global discos;
parameter_optimized_Decision_Surface=true;
use_accuracy=true;
Currentfolder=pwd;
addpath(genpath(fullfile(Currentfolder,'Classifiers','SLR','Classifier_Files')));

%Load some paramters about the RCE algorithm from the Results file.
global discos folder_name;
cd(fullfile(folder_name,'Rce_Results',discos));
load Results.mat
no_iter =kfold*topResample;
%clearvars -except  clas9 CurrentFolder test_dat kfold topResample classes classname no_iter RCE_steps parameter_optimized_Decision_Surface max_accuracy_Decision_Surface use_accuracy acc_for_each_iteration acc_for_each_iteration_opt FinalConfmat FinalConfmat_opt 
cd(Currentfolder);
Currentfolder=pwd;
% Load test data
%load test_data.mat 
%[~,~,test_data]=xlsread('train.xlsx');
%test_data=test_dat;
test_data=handles.test_data;
test_data_targets = test_data(2,3:end); 
coordinates = test_data(3:end,2);
test_data=test_data(3:end,3:end);
test_data_inp =(cell2mat(test_data))';

% Convert group names to group nos for the testing data 
no_classes=size(classname,2);
group_no_targets=zeros(size(test_data_inp,1),1);
for ns=1:1:size(test_data_inp,1);
    for ng=1:no_classes
        if(strcmp(test_data_targets(ns),classname{ng})==1)
        group_no_targets(ns)=ng;
        end
    end
end
cd(fullfile(folder_name,'Rce_Results',discos));
% Predict the values of  Test data

% Concatenate the predictions from all decision Surfaces
cd Decision_Surfaces
h = waitbar(0,'Initializing waitbar...');
if ~max_accuracy_Decision_Surface
    Confmat_test=zeros(no_classes,no_classes,RCE_steps);
    list = dir('Decision_Surface*.mat');
    classes_predicted =[]; test_accuracy_ite =[];
    for k=1:length(list)
        load(list(k).name)  
        class_act = zeros(size(Decision_surface,1),size(Decision_surface,2),size(test_data_inp,1));
        test_accuracy_resample =zeros(size(Decision_surface,1),size(Decision_surface,2));
     for i=1:size(Decision_surface,2)       
         for j=1:size(Decision_surface,1)
             if parameter_optimized_Decision_Surface
                  model = Decision_surface_optimized{j,i};
             else
                  model=Decision_surface{j,i};
             end
              [class_act(j,i,:)] = classifier_predict(model,test_data_inp(:,Spot_ID{j,i}));
              test_accuracy_resample(j,i) =mean(squeeze(class_act(j,i,:))==group_no_targets);
         end 
     end
     waitbar(k/topResample,h,sprintf('Calculating... %d%%',k*100/topResample));
          clear Decision_Surface Spot_ID_TF
         classes_predicted = cat(1,classes_predicted, class_act);
         test_accuracy_ite = cat(1,test_accuracy_ite,test_accuracy_resample); 
    end
    test_data_accuracy_voting = zeros(1, RCE_steps); test_data_specificity_voting = zeros(1, RCE_steps);
    test_data_sensitivity_voting = zeros(1, RCE_steps); 
    for r =1: RCE_steps
         cp = classperf(group_no_targets);
         test_data_predict = mode(squeeze(classes_predicted(:,r,:)),1);
         classperf(cp, test_data_predict);
         test_data_accuracy_voting(r)= cp.CorrectRate;
         test_data_specificity_voting(r)= cp.Specificity;
         test_data_sensitivity_voting(r)=cp.Sensitivity;
         Confmat_test(:,:,r)=confusionmat( group_no_targets,test_data_predict,'order',unique(group_no_targets));
    end
         test_data_accuracy_actual = mean(test_accuracy_ite,1);
         test_data_accuracy_std = std(test_accuracy_ite,1);
        if parameter_optimized_Decision_Surface
            CV_accuracy = mean(acc_for_each_iteration_opt,1);
            CV_accuracy_std = std(acc_for_each_iteration_opt,1);
            Confmat_CV = FinalConfmat_opt;
        else
            CV_accuracy = mean(acc_for_each_iteration,1);
            CV_accuracy_std = std(acc_for_each_iteration,1);
            Confmat_CV = FinalConfmat;
        end
       
        for j=1:no_classes
             assignin('base',strcat('acc_CV_',classname{j}),(squeeze(Confmat_CV(j,j,:)))'./sum(squeeze(Confmat_CV(j,:,:)),1));
             assignin('base',strcat('acc_test_',classname{j}),(squeeze(Confmat_test(j,j,:)))'./sum(squeeze(Confmat_test(j,:,:)),1));
        end
        
else     % If only the RCE step with the best performance is chosen
 if use_accuracy && parameter_optimized_Decision_Surface
     cd Accuracy_optimized
 elseif use_accuracy
     cd Accuracy;
 else
    cd Performance;
 end
    list = dir('Decision_Surface*.mat');
    cp = classperf(group_no_targets);
    classes_predicted=[]; test_accuracy_ite =[];
    for k=1:length(list)
        load(list(k).name)  
        test_accuracy_resample = zeros(size(Decision_surface,1),1);
        class_act = zeros(size(Decision_surface,1),size(test_data_inp,1));      
        for j=1:size(Decision_surface,1)
             model = Decision_surface{j};
             index = ismember(coordinates,Paths{j});
             test_data_red = test_data_inp(:,index);
             coordinates_red = coordinates(index);
             [class_act(j,:)] = classifier_predict(model,test_data_red);
             test_accuracy_resample(j) =mean(class_act(j,:)== group_no_targets');
        end
         waitbar(k/topResample,h,sprintf('Calculating... %d%%',k*100/topResample));
         clear Decision_Surface* index coordinates_red
         classes_predicted = cat(1,classes_predicted, class_act);
         test_accuracy_ite = cat(1,test_accuracy_ite,test_accuracy_resample); 
    end
         test_data_predict = mode(classes_predicted,1);
         classperf(cp,test_data_predict);
         test_data_accuracy_actual = mean(test_accuracy_ite);
         test_data_accuracy_std = std(test_accuracy_ite);
         test_data_accuracy_voting= cp.CorrectRate;
         test_data_specificity_voting = cp.Specificity;
         test_data_sensitivity_voting = cp.Sensitivity;
         Confmat_test = confusionmat(group_no_targets,test_data_predict,'order',unique(group_no_targets));
         disp(Confmat_test);
         cd ..
        clear Decision_surface*
        if parameter_optimized_Decision_Surface
            [CV_accuracy, Indx] = max(mean(acc_for_each_iteration_opt,1));
            CV_accuracy_std = std(squeeze(acc_for_each_iteration_opt(:,Indx)),1);
            Confmat_CV = squeeze(FinalConfmat_opt(:,:,Indx));
        else
            [CV_accuracy, Indx] = max(mean(acc_for_each_iteration,1));
            CV_accuracy_std = std(squeeze(acc_for_each_iteration(:,Indx)));
            Confmat_CV = squeeze(FinalConfmat(:,:,Indx));
        end
        for j=1:no_classes
            temp(2*(j-1)+1) = Confmat_CV(j,j)/sum(squeeze(Confmat_CV(j,:)));
            temp(2*j) = Confmat_test(j,j)/sum(squeeze(Confmat_test(j,:)));
            assignin('base',strcat('acc_CV_',classname{j}),Confmat_CV(j,j)/sum(squeeze(Confmat_CV(j,:))));
            assignin('base',strcat('acc_test_',classname{j}), Confmat_test(j,j)/sum(squeeze(Confmat_test(j,:))));
        end
                end
balanced_CV_accuracy=mean(temp(1:2:end));
balanced_test_accuracy=mean(temp(2:2:end));
rmpath(genpath(fullfile(Currentfolder,'Classifiers','SLR','Classifier_Files')));
cd(fullfile(folder_name,'Rce_Results',discos));
close(h);
save Test_results -regexp ^(?!(CurrentFolder|test_dat|kfold|topResample|classes|classname|no_iter|RCE_steps|parameter_optimized_Decision_Surface|max_accuracy_Decision_Surface|use_accuracy|acc_for_each_iteration|acc_for_each_iteration_opt|FinalConfmat|FinalConfmat_opt|time|timer|tmse|TOP_RESAMPLE|total_freq|ts_fmri_f|TSfiles|uppe|v|X1|X2|y1|y2|YourFile|clas2|clas3|clas4|clas5|clas6|clas7|clas8|clas9|clas10|clas11|clas12|clas13|clas14|clas15|clas16|clas17|clas18|clas19|clas20|clas21|classes_predicted|data|Data|data_test|data_train|data_train_tmp|Decision_name|Decision_surf_ite|Decision_surf_opt_ite|decrease|desi|desicoll|dfg|DSF|end_nmc|end_paths|Estimate_size_time|eventdata|fid|fileID|Filepath|first|folder_name|fout|fpaths|freq|fsamples|gfd|gindex|group_no_targets|Groupdir|groups_ID|groups_name|groups_test|groups_train|grp|gsorted|h|handles|hObject|hps|HPS|Hyperparameter|Hyperparameter_Search|i|idod|iip|IIP|incred|indices2|indx|init_nmc|input_data|j|jj|k|kk|l|Labels|last_num_clusters|list|loc|lowv|madc|MADC|lowto|incredto|uppeto|n|n2|n3|n_samples|ng|nmc|no_classes|no_clusters|nocl|nofo|ns|resa|resample|row|sadi|score|sdir|select_significant_paths|selected_paths_resample|significant_paths|specificity|specificity_opt|splirat|split_point|Spot_id|Spot_ID|Spot_ID_acc|Spot_ID_acc_opt|Spot_ID_ite|Spot_ID_perf|Spot_ID_TF|Spot_ID_tmp|Spotid|ssp|SSP|st|staticFC_Controls|staticFC_PTSD|stfc|sth|Table|temp|test)$).
cd(Currentfolder)
end