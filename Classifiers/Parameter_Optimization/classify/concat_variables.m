function [accuracy_final,specificity_final, sensitivity_final, ConfusionMatrix_final, varargout] = concat_variables(CurrentFolder,topResample,kfold,no_classes,num_paths,no_Hyperparameters,HyperParameter_Search,identify_important_paths)
global discos folder_name;
% initialize variables
ConfusionMatrix_final = zeros(no_classes,no_classes,kfold*topResample); accuracy_final = zeros(kfold*topResample,1); 
specificity_final= zeros(kfold*topResample,1); sensitivity_final = zeros(kfold*topResample,1);
variable_importance_score_final = zeros(num_paths,kfold*topResample); 
if HyperParameter_Search
    best_parameters_final = zeros(kfold*topResample,no_Hyperparameters);
end
mkdir(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Accuracy'));
cd(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Accuracy'));
listing = dir(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Accuracy',filesep,'*.mat'));
 for list=1:topResample
     load(listing(list).name)
     accuracy_final((kfold*(list-1)+1):(kfold*list))  = accuracy_allclass;
     specificity_final((kfold*(list-1)+1):(kfold*list)) = specificity_allclass;
     sensitivity_final((kfold*(list-1)+1):(kfold*list)) = sensitivity_allclass;
     ConfusionMatrix_final(:,:,(kfold*(list-1)+1):(kfold*list)) = Conf_mattest;
 end
 clear listing 
   if HyperParameter_Search
        mkdir(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'bestparameters'));
       cd(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'bestparameters'));
       listing = dir(fullfile(folder_name, filesep,'Outside_Rce_Results',discos,filesep,'bestparameters',filesep,'*.mat'));
       for list = 1:topResample
           load(listing(list).name)
           best_parameters_final((kfold*(list-1)+1):(kfold*list),:)  = best_parameters;
       end
    clear listing best_parameters;
    cd(fullfile(folder_name,filesep,'Outside_Rce_Results',discos))
    save('bestparameters.mat','best_parameters_final')
   end

if identify_important_paths
    mkdir(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Feature_Significance'));
     cd(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Feature_Significance'));
   listing = dir(fullfile(folder_name, filesep,'Outside_Rce_Results',discos,filesep,'Feature_Significance',filesep,'*.mat'));
   for list = 1:topResample
           load(listing(list).name)
           variable_importance_score_final(:,(kfold*(list-1)+1):(kfold*list)) = variable_importance_score;
   end
   clear listing variable_importance_score;
   varargout{1}=sum(variable_importance_score_final,2);
end

