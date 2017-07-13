function varargout = Run_Classifier(CurrentFolder,input_data,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,HyperParameter_Search,identify_important_paths,topResample,Estimate_size_time,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input_data = ReadFromExcel(filename,'ALL');
global discos folder_name;
if nargin > 11
  Hyperparameter =varargin{1};
else
  Hyperparameter={};
end

[row, col] = size(input_data);         
data = cell2mat(input_data(3:row,3:col));
samples_name  = cell2mat(input_data(1,3:col));                % Extract subjects 
[~, n_samples]= size(samples_name);
groups_name   = input_data(2,3:col);                % Extract the group names
paths = input_data(3:row,2);
spot_id = cell2mat(input_data(3:row,1));
classes = 1:length(unique(groups_name));
% Convert Group names into Group Indices
groups_ID=zeros(n_samples,1);
for ns=1:n_samples
    for ng=1:no_classes
        if(strcmp(groups_name(ns),classname{ng})==1)
        groups_ID(ns)=ng;
        end
    end
end

if size(unique(groups_ID),1) ~= no_classes
    error('The no of unique groups in the data do not match. Please check');
end

clear input_data;                                     %free memory

[num_paths, num_subjects] = size(data); 


timer =tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run k-folds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data= data';                                      %treat the samples in rows
TOP_RESAMPLE = topResample;
outer_kfold = kfold;
for resample =1:TOP_RESAMPLE
    fprintf('\n**************Top Level Resample=%d*************',resample);
    indices2 = crossvalind('Kfold',groups_ID,outer_kfold);
    Decision_surface = cell(outer_kfold,1);
    accuracy_allclass =zeros(outer_kfold,1);
    specificity_allclass = zeros(outer_kfold,1);
    sensitivity_allclass = zeros(outer_kfold,1);
    Conf_mattest = zeros(no_classes, no_classes,outer_kfold);
     if select_significant_paths
        selected_paths_resample = zeros(num_paths,outer_kfold);
          Paths_selected= cell(outer_kfold,1);
     end
if HyperParameter_Search
    best_parameters = zeros(outer_kfold,length(Hyperparameter));
end
if identify_important_paths
    variable_importance_score = zeros(num_paths,outer_kfold);
end 
for i=1:outer_kfold
    fprintf('\n#Iteration:%d\n',((resample-1)*outer_kfold +i));
 
     spotid=spot_id;                                %The data is represented by paths at columns and samples at row
     test  = (indices2 == i); 
     train =~ test;
     data_train   = data(train,:);
     data_test    = data(test,:);
     groups_train = groups_ID(train);
     groups_test  = groups_ID(test);
     
%  Multivatiate N-way Anova (i.e. N-way hypothesis testing) 

       data_train_tmp = data_train;
       data_test_tmp  = data_test;
       paths_tmp = paths;  
       spot_id_tmp=spot_id ; 
   if select_significant_paths
       [mask] = select_significant_paths_fun(data_train_tmp,classes,groups_train);
       data_train_tmp = data_train_tmp(:,mask);
       data_test_tmp  = data_test_tmp(:,mask);
       Paths_ite   = paths_tmp(mask);
       spot_id_tmp = spot_id_tmp(mask);
       Paths_selected {i} =Paths_ite;
    end
  
%       data_train_tmp = zscore(data_train_tmp')';
%       data_test_tmp  = zscore (data_test_tmp')';

     selected_paths_resample(:,i) = ismember(spotid,spot_id_tmp);     
 if HyperParameter_Search

    [accuracy_allclass(i), specificity_allclass(i), sensitivity_allclass(i),Conf_mattest(:,:,i), Decision_surf_ite,best_parameters(i,:)] = optimPara(data_train_tmp,groups_train, data_test_tmp,groups_test,kfold,1,classes,Hyperparameter);

 else 
     
    [accuracy_allclass(i), specificity_allclass(i), sensitivity_allclass(i), Conf_mattest(:,:,i), Decision_surf_ite] = runclassifier_Default( data_train_tmp,groups_train,data_test_tmp,groups_test,classes);

 end

 
 disp(squeeze(Conf_mattest(:,:,i)))
 if identify_important_paths
     variable_importance_score(:,i) = find_predictorImportance(Decision_surf_ite,spot_id_tmp,spotid);
 end
 
  if Decision_Surface_Flag
     Decision_surface{i} = Decision_surf_ite;
  end 

end
mkdir(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Accuracy'));
cd(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Accuracy'));
savevars(sprintf('accuracyvals%d.mat', resample), 'accuracy_allclass',accuracy_allclass, 'specificity_allclass',specificity_allclass,'sensitivity_allclass',sensitivity_allclass,'Conf_mattest',Conf_mattest);

 if Decision_Surface_Flag
     mkdir(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Decision_Surfaces')); 
    cd(fullfile(folder_name,filesep,'Outside_Rce_Results',discos,filesep,'Decision_Surfaces'));  
    Decision_name = strcat('Decision_Surface',num2str(resample),'.mat');
    if select_significant_paths
         savevars(Decision_name,'Decision_Surface',Decision_surface, 'paths_selected',selected_paths_resample,'Paths', Paths_selected);
    else
        savevars(Decision_name,'Decision_Surface',Decision_surface);
    end
 end
 
 if HyperParameter_Search
    mkdir(fullfile(folder_name, filesep,'Outside_Rce_Results',discos,filesep,'bestparameters'))
     cd(fullfile(folder_name, filesep,'Outside_Rce_Results',discos,filesep,'bestparameters'))
     savevars(sprintf('bestparams%d.mat', resample),'best_parameters',best_parameters);
 end
 if identify_important_paths
     mkdir(fullfile(folder_name, filesep,'Outside_Rce_Results',discos,filesep,'Feature_Significance'))
     cd(fullfile(folder_name, filesep,'Outside_Rce_Results',discos,filesep,'Feature_Significance'))
     savevars(sprintf('variable_importance%d.mat', resample),'variable_importance_score',variable_importance_score);
 end
 
 cd(fullfile(CurrentFolder,'Classifiers','Parameter_Optimization','classify'))
end
if ~Estimate_size_time
if identify_important_paths
  [accuracy_final,specificity_final, sensitivity_final, ConfusionMatrix_final, variable_importance_sum] = concat_variables(CurrentFolder, TOP_RESAMPLE,outer_kfold,no_classes,num_paths,length(Hyperparameter),HyperParameter_Search,identify_important_paths);
  [variable_importance_vals, variable_importance_order] = sort(variable_importance_sum,'descend');
else  
  [accuracy_final, specificity_final, sensitivity_final, ConfusionMatrix_final] = concat_variables(CurrentFolder, TOP_RESAMPLE,outer_kfold,no_classes,num_paths,length(Hyperparameter),HyperParameter_Search,identify_important_paths);
end


FinalConfmat =sum(ConfusionMatrix_final,3);
[precision,recall,specificity]=acc_cal(FinalConfmat);

cd(fullfile(folder_name,'Outside_Rce_Results',discos));
mean_acc = mean (accuracy_final);
stdv_acc = std (accuracy_final);

fout =  fopen('Final Accuracy Results.txt','w');
fprintf(fout, 'Final Accuracy %3.2f\n',mean_acc); 
fprintf(fout, 'Std Dev Accuracy %3.2f\n\n',stdv_acc); 
accuracy = zeros(no_classes,1);
for j=1:no_classes
accuracy(j) = (FinalConfmat(j,j)/sum(squeeze(FinalConfmat(j,:))));
end
for i=1:no_classes
 fprintf(fout,'Class Name = %s\n',classname{i});
 fprintf(fout,'Accuracy %3.2f \tPrecision %3.2f \tRecall %3.2f \tSpecificity %3.2f\n\n',accuracy(i),precision(i),recall(i),specificity(i));
end
fclose(fout);
end
if Decision_Surface_Flag && Estimate_size_time
 cd(fullfile(folder_name,'Outside_Rce_Results',discos,'Decision_Surfaces'));
 S= dir('Decision_Surface1.mat');
 Decision_Surf_size = S.bytes;  varargout{2} = Decision_Surf_size;
end
cd(fullfile(CurrentFolder,'Classifiers','Parameter_Optimization'))
if ~Estimate_size_time
    cd ..
    cd(fullfile(folder_name,'Outside_Rce_Results',discos))
    save('Results.mat', '-v7.3'); %will save the workspace into Results_sfc.mat
end
time=toc(timer); varargout{1}=time;
cd(fullfile(CurrentFolder,'Classifiers','Parameter_Optimization','classify'))
 