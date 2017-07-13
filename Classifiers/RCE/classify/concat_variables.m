function [FinalConfmat,cum_perf_Final,acc_for_each_iteration,max_performance_per_cluster,performance_each_cluster_test,performance_each_cluster_train,...
num_paths_for_max_clusters,num_paths_iter,num_cluster_iter,all_paths_freq,varargout] = concat_variables(CurrentFolder,HyperParameter_Search, max_accuracy_Decision_Surface,no_classes,RCE_steps,kfold,topResample,paths,varargin)
global discos folder_name;    

if   HyperParameter_Search
        no_Hyperparameters =varargin{1};
    end
        mkdir(fullfile(folder_name,'Rce_Results',discos,'Accuracy_vals'))
   cd(fullfile(folder_name,'Rce_Results',discos,'Accuracy_vals'))
   listing = dir(fullfile(folder_name,'Rce_Results',discos,'Accuracy_vals','*.mat'));
   num_paths_iter =zeros (kfold*topResample,RCE_steps); FinalConfmat = zeros(no_classes,no_classes,RCE_steps); cum_perf_Final = zeros(4, no_classes, RCE_steps);
   acc_for_each_iteration = zeros(kfold*topResample, RCE_steps); max_performance_per_cluster = zeros(kfold*topResample, RCE_steps);
   performance_each_cluster_test = zeros(kfold*topResample, RCE_steps); performance_each_cluster_train =zeros(kfold*topResample, RCE_steps);
   num_paths_for_max_clusters = zeros(kfold*topResample, RCE_steps); num_cluster_iter =zeros(kfold*topResample, RCE_steps);
   all_paths_freq  =  zeros(length(paths),RCE_steps);  max_accuracy_all_iter = zeros(kfold*topResample,1); max_performance_all_iter = zeros(kfold*topResample,1); 
   if HyperParameter_Search
       FinalConfmat_opt = zeros(no_classes,no_classes,RCE_steps); cum_perf_Final_opt=zeros(4, no_classes, RCE_steps);acc_for_each_iteration_opt=zeros(kfold*topResample, RCE_steps);
       best_parameters_final = zeros(kfold*topResample,no_Hyperparameters);
   end    
   for list=1:topResample
       load(listing(list).name)
       FinalConfmat = FinalConfmat+Confmat_resample;
       cum_perf_Final =  cum_perf_Final+cum_perf_resample;
       acc_for_each_iteration((kfold*(list-1)+1):(kfold*list),:) = acc_for_each_resample;
       max_performance_per_cluster((kfold*(list-1)+1):(kfold*list),:) = max_performance_per_cluster_resample;
       performance_each_cluster_test((kfold*(list-1)+1):(kfold*list),:) = performance_each_resample_test;
       performance_each_cluster_train((kfold*(list-1)+1):(kfold*list),:) = performance_each_resample_train;
       num_paths_for_max_clusters((kfold*(list-1)+1):(kfold*list),:) = num_paths_for_best_clusters_resample;
       num_paths_iter((kfold*(list-1)+1):(kfold*list),:) = num_paths_resample;
       num_cluster_iter((kfold*(list-1)+1):(kfold*list),:) = last_num_clusters;
       if  HyperParameter_Search
           FinalConfmat_opt = FinalConfmat_opt + Confmat_resample_opt;
           cum_perf_Final_opt =  cum_perf_Final_opt+cum_perf_resample_opt;
           acc_for_each_iteration_opt((kfold*(list-1)+1):(kfold*list),:) = acc_for_each_resample_opt;
           best_parameters_final((kfold*(list-1)+1):(kfold*list),:)  = best_parameters;
       end
       
       for i =1:kfold 
         paths_freq_fold = paths_freq_resample{i};
         selected_paths =selected_paths_resample(:,i);
        indx=0; 
        for q =1:size(selected_paths,1)
            if (selected_paths(q)==1)
                indx = indx + 1;
                all_paths_freq(q,:) = all_paths_freq(q,:) + paths_freq_fold(indx,:); %Keep the frequency of each path
            end
        end
       end
   end
 if max_accuracy_Decision_Surface
      mkdir(fullfile(folder_name,'Rce_Results',discos,'MaxAccuracy_vals'))
  cd(fullfile(folder_name,'Rce_Results',discos,'MaxAccuracy_vals'))
  listing = dir(fullfile(folder_name,'Rce_Results',discos,'MaxAccuracy_vals','*.mat'));
   for list=1:length(listing)
       load(listing(list).name)
       max_accuracy_all_iter((kfold*(list-1)+1):(kfold*list)) = max_accuracy_each_iter;
       max_performance_all_iter((kfold*(list-1)+1):(kfold*list)) = max_performance_each_iter;
       delete(listing(list).name)
   end
  cd ..
  save('max_Accuracy.mat','max_accuracy_all_iter', 'max_performance_all_iter')
  %rmdir(fullfile(CurrentFolder,'Results','MaxAccuracy_vals'))
 end
  if HyperParameter_Search
     cd(fullfile(folder_name,filesep,'Rce_Results',discos))
    save('bestparameters.mat','best_parameters_final')
  end

  if HyperParameter_Search
 varargout{1}= FinalConfmat_opt; varargout{2}= cum_perf_Final_opt; varargout{3}=acc_for_each_iteration_opt;
  end