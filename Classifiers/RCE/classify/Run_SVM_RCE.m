
function varargout = Run_SVM_RCE(CurrentFolder,input_data,init_nmc,end_nmc,decrease,end_paths,kfold,classname,no_classes,select_significant_paths,Decision_Surface_Flag,max_accuracy_Decision_Surface,HyperParameter_Search,topResample,Estimate_size_time, varargin)
global discos folder_name;

% Parameters:
% init_nmc : the number of initialized clusters
% end_nmc  : the number of clusters at the final stage
% decrease : decrease in each stage (type in integer,e.g. 0.1=10%)
% end_paths : the number of paths at the final stage
% topResample: the maximum number of the resamples
% classname: array containing the target class names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input file prototype:
%   SpotID Pathname  Sample_1 Sample_2 ..... Sample_N
%    ID1   Name1     Class_1  Class_1  ..... Class_2
%    ID2   Name2     values   values   ..... values
%     .      .         .         .       .      .       
%     .      .         .         .       .      .
%     .      .         .         .       .      .
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input_data = ReadFromExcel(filename,'ALL');
if nargin > 15
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


%find the number of elements in each classes

timer=tic;

%find the number of RCE Steps
RCE_steps = 0;
nmc=init_nmc;
while nmc >= end_nmc 
RCE_steps =RCE_steps+1; 
no_clusters(RCE_steps)=nmc;
nmc = floor(nmc*(1-decrease));
end

% initialize variables
if Decision_Surface_Flag && max_accuracy_Decision_Surface && Estimate_size_time
    mkdir(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces','Accuracy'))
    mkdir(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces','Performance'))
    if HyperParameter_Search
       mkdir(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces','Accuracy_optimized'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run k-folds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = data';                                      %treat the samples in rows
TOP_RESAMPLE = topResample;
for resample =1:TOP_RESAMPLE
    performance_each_resample_test = zeros(kfold,RCE_steps); performance_each_resample_train = zeros(kfold,RCE_steps); max_performance_per_cluster_resample = zeros(kfold,RCE_steps);
    num_paths_for_best_clusters_resample = zeros(kfold,RCE_steps); cum_perf_resample = zeros(4, no_classes, RCE_steps); Confmat_resample = zeros(no_classes,no_classes,RCE_steps);
    acc_for_each_resample =zeros(kfold,RCE_steps); max_accuracy_each_iter=[]; max_performance_each_iter=[]; paths_freq_resample= []; num_paths_resample = zeros(kfold, RCE_steps);
    last_num_clusters = zeros(kfold, RCE_steps); selected_paths_resample = zeros(num_paths, kfold);
    if HyperParameter_Search 
        acc_for_each_resample_opt =zeros(kfold,RCE_steps); cum_perf_resample_opt = zeros(4, no_classes, RCE_steps); Confmat_resample_opt =zeros(no_classes,no_classes,RCE_steps);
        max_accuracy_each_iter_opt=[]; best_parameters = zeros(kfold,length(Hyperparameter));
    end   
    if Decision_Surface_Flag
        num_paths_TF=[]; Paths_TF = []; Spot_ID_TF =[]; Decision_surface=[];
        if HyperParameter_Search 
           Decision_surface_opt=[]; 
        end
   end
    fprintf('\n**************Top Level Resample=%d*************',resample);
    indices2 = crossvalind('Kfold',groups_ID,kfold);
    outer_kfold = kfold;  
    for i=1:outer_kfold        
        fprintf('\n#Iteration:%d\n',((resample-1)*outer_kfold +i));
        paths_tmp = paths;  
        spot_id_tmp=spot_id ;   
        spotid=spot_id;                                %The data is represented by paths at columns and samples at row
         test  = (indices2 == i); 
         train =~ test;
         data_train   = data(train,:);
         data_test    = data(test,:);
         groups_train = groups_ID(train);
         groups_test  = groups_ID(test);
         samples_train= samples_name(train);
         samples_test = samples_name(test);
         
 %  Multivatiate N-way Anova (i.e. N-way hypothesis testing)          
       data_train_tmp = data_train;
       data_test_tmp  = data_test;
       paths_tmp = paths;  
       spot_id_tmp=spot_id ; 
       
   if select_significant_paths 
       [mask] = select_significant_paths_fun(data_train_tmp,classname,groups_train);
       data_train_tmp = data_train_tmp(:,mask);
       data_test_tmp  = data_test_tmp(:,mask);
       paths_tmp   = paths_tmp(mask);
       spot_id_tmp = spot_id_tmp(mask);
    end
  
%       data_train_tmp = zscore(data_train_tmp')';
%       data_test_tmp  = zscore (data_test_tmp')';

          selected_paths_resample(:,i) = ismember(spotid,spot_id_tmp);                     %Point to the paths been selected to the analysis
        
         if HyperParameter_Search && Decision_Surface_Flag
            
            [max_avg_performance_each_cluster, Conf_mattest, performance, last_num_clusters(i,:), num_paths_resample(i,:), paths_freq, Decision_surf_ite, num_paths_ite, Paths_ite, Spot_ID_ite, performance_opt, Conf_mattest_opt, Decision_surf_opt_ite,best_parameters(i,:)] = core_svmKmeansGeneSelection_v3( data_train_tmp',spot_id_tmp,paths_tmp,groups_train,samples_train,...
            init_nmc,end_nmc,end_paths,decrease,kfold,1,...
            classes,data_test_tmp',groups_test,test,RCE_steps,Decision_Surface_Flag,HyperParameter_Search,Hyperparameter);
             num_paths_TF=cat(1,num_paths_TF,num_paths_ite); Paths_TF = cat(1,Paths_TF,Paths_ite); Spot_ID_TF =cat(1,Spot_ID_TF,Spot_ID_ite); Decision_surface=cat(1,Decision_surface,Decision_surf_ite); Decision_surface_opt=cat(1,Decision_surface_opt,Decision_surf_opt_ite);
        elseif HyperParameter_Search && ~Decision_Surface_Flag 
            
            [max_avg_performance_each_cluster, Conf_mattest, performance, last_num_clusters(i,:), num_paths_resample(i,:), paths_freq, performance_opt, Conf_mattest_opt,best_parameters(i,:)] = core_svmKmeansGeneSelection_v3( data_train_tmp',spot_id_tmp,paths_tmp,groups_train,samples_train,...
            init_nmc,end_nmc,end_paths,decrease,kfold,1,...
            classes,data_test_tmp',groups_test,test,RCE_steps,Decision_Surface_Flag,HyperParameter_Search,Hyperparameter);

        elseif ~HyperParameter_Search && Decision_Surface_Flag
            
            [max_avg_performance_each_cluster, Conf_mattest, performance, last_num_clusters(i,:), num_paths_resample(i,:), paths_freq, Decision_surf_ite, num_paths_ite, Paths_ite, Spot_ID_ite] = core_svmKmeansGeneSelection_v3( data_train_tmp',spot_id_tmp,paths_tmp,groups_train,samples_train,...
            init_nmc,end_nmc,end_paths,decrease,kfold,1,...
            classes,data_test_tmp',groups_test,test,RCE_steps,Decision_Surface_Flag,HyperParameter_Search);
            num_paths_TF=cat(1,num_paths_TF,num_paths_ite); Paths_TF = cat(1,Paths_TF,Paths_ite); Spot_ID_TF =cat(1,Spot_ID_TF,Spot_ID_ite); Decision_surface=cat(1,Decision_surface,Decision_surf_ite);
        else       
            
            [max_avg_performance_each_cluster, Conf_mattest, performance, last_num_clusters(i,:), num_paths_resample(i,:), paths_freq] = core_svmKmeansGeneSelection_v3( data_train_tmp',spot_id_tmp,paths_tmp,groups_train,samples_train,...
            init_nmc,end_nmc,end_paths,decrease,kfold,1,...
            classes,data_test_tmp',groups_test,test,RCE_steps,Decision_Surface_Flag,HyperParameter_Search);

         end
         paths_freq_resample{i} = paths_freq;
         
        max_performance_per_cluster_resample(i,:) = squeeze(max_avg_performance_each_cluster(3,:)); %cum max performance on cluster on training data
        performance_each_resample_test(i,:) = squeeze(max_avg_performance_each_cluster(1,:));
        performance_each_resample_train(i,:) = squeeze(max_avg_performance_each_cluster(2,:));       
        num_paths_for_best_clusters_resample(i,:)  = max_avg_performance_each_cluster(4,:); %Avg of paths at the maximum clusters
        cum_perf_resample = cum_perf_resample+ performance;
        Confmat_resample= Confmat_resample+ Conf_mattest;
        
        one_itr_acc= sum(performance(1,:,:),2)./(sum(sum(Conf_mattest,1),2)); % accuracy per iteration
        acc_for_each_resample(i,:) = squeeze(one_itr_acc);
        
        if HyperParameter_Search   
        cum_perf_resample_opt =  cum_perf_resample_opt+performance_opt;
        Confmat_resample_opt = Confmat_resample_opt + Conf_mattest_opt;
        one_itr_acc_opt= sum(performance_opt(1,:,:),2)./(sum(sum(Conf_mattest_opt,1),2)); % accuracy per iteration
        acc_for_each_resample_opt(i,:) = squeeze(one_itr_acc_opt);
        end
    end
        mkdir(fullfile(folder_name,'Rce_Results',discos,'Accuracy_vals')) 
      cd(fullfile(folder_name,'Rce_Results',discos,'Accuracy_vals')) 
      if HyperParameter_Search
          savevars(sprintf('AccuracyResults%d.mat', resample),'max_performance_per_cluster_resample', max_performance_per_cluster_resample,'performance_each_resample_test',performance_each_resample_test,...
         'performance_each_resample_train', performance_each_resample_train,'num_paths_for_best_clusters_resample', num_paths_for_best_clusters_resample,'cum_perf_resample', cum_perf_resample,'Confmat_resample',...
         Confmat_resample,'acc_for_each_resample',acc_for_each_resample,'selected_paths_resample',selected_paths_resample,'last_num_clusters',last_num_clusters, 'num_paths_resample',num_paths_resample, 'paths_freq_resample', paths_freq_resample,'cum_perf_resample_opt',cum_perf_resample_opt,...
         'Confmat_resample_opt', Confmat_resample_opt,'acc_for_each_resample_opt',acc_for_each_resample_opt,'best_parameters',best_parameters);
      else
         savevars(sprintf('AccuracyResults%d.mat', resample),'max_performance_per_cluster_resample', max_performance_per_cluster_resample,'performance_each_resample_test',performance_each_resample_test,...
         'performance_each_resample_train', performance_each_resample_train,'num_paths_for_best_clusters_resample', num_paths_for_best_clusters_resample,'cum_perf_resample', cum_perf_resample,'Confmat_resample',...
         Confmat_resample,'acc_for_each_resample',acc_for_each_resample,'selected_paths_resample',selected_paths_resample,'last_num_clusters',last_num_clusters, 'num_paths_resample',num_paths_resample, 'paths_freq_resample', paths_freq_resample);
     end
     
      if Decision_Surface_Flag
          mkdir(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces'));
       cd(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces'));
       Decision_name = strcat('Decision_Surface',num2str(resample),'.mat');
       if max_accuracy_Decision_Surface
             if HyperParameter_Search
                [accuracy_steps, RCE_steps_sort] = sort(acc_for_each_resample_opt,2,'ascend');
                max_accuracy_opt= accuracy_steps(:,end); RCE_step_opt =RCE_steps_sort(:,end);
                num_paths_acc_opt  = (num_paths_TF(sub2ind(size(num_paths_TF),1: size(num_paths_TF,1),RCE_step_opt')))';
                Paths_acc_opt = (Paths_TF(sub2ind(size(Paths_TF),1: size(Paths_TF,1),RCE_step_opt')))';
                Spot_ID_acc_opt = (Spot_ID_TF(sub2ind(size(Spot_ID_TF),1: size(Spot_ID_TF,1),RCE_step_opt')))';
                Decision_surface_acc_opt =(Decision_surface_opt(sub2ind(size(Decision_surface_opt),1: size( Decision_surface_opt,1),RCE_step_opt')))';
                max_accuracy_each_iter_opt = cat(1,max_accuracy_each_iter_opt,max_accuracy_opt);
             end
             [accuracy_steps, RCE_steps_sort] = sort(acc_for_each_resample,2,'ascend');
             max_accuracy= accuracy_steps(:,end); RCE_step =RCE_steps_sort(:,end);
             num_paths_acc  = (num_paths_TF(sub2ind(size(num_paths_TF),1: size(num_paths_TF,1),RCE_step')))';
             Paths_acc = (Paths_TF(sub2ind(size(Paths_TF),1: size(Paths_TF,1),RCE_step')))';
             Spot_ID_acc = (Spot_ID_TF(sub2ind(size(Spot_ID_TF),1: size(Spot_ID_TF,1),RCE_step')))';
             Decision_surface_acc =(Decision_surface(sub2ind(size(Decision_surface),1: size( Decision_surface,1),RCE_step')))';
             max_accuracy_each_iter = cat(1,max_accuracy_each_iter,max_accuracy);     

             [accuracy_steps, RCE_steps_sort] = sort(performance_each_resample_test,2,'ascend');
             max_performance = accuracy_steps(:,end); RCE_step =RCE_steps_sort(:,end);
             num_paths_perf  = (num_paths_TF(sub2ind(size(num_paths_TF),1: size(num_paths_TF,1),RCE_step')))';
             Paths_perf = (Paths_TF(sub2ind(size(Paths_TF),1: size(Paths_TF,1),RCE_step')))';
             Spot_ID_perf = (Spot_ID_TF(sub2ind(size(Spot_ID_TF),1: size(Spot_ID_TF,1),RCE_step')))';
             Decision_surface_perf =(Decision_surface(sub2ind(size(Decision_surface),1: size( Decision_surface,1),RCE_step')))';      
             max_performance_each_iter = cat(1,max_performance_each_iter,max_performance); 
             mkdir(fullfile(folder_name,'Rce_Results',discos,'MaxAccuracy_vals'));
             cd(fullfile(folder_name,'Rce_Results',discos,'MaxAccuracy_vals'));
             savevars(sprintf('max_Accuracy%d.mat', resample),'max_performance_each_iter',max_performance_each_iter,'max_accuracy_each_iter',max_accuracy_each_iter);
             if HyperParameter_Search
                 mkdir(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces','Accuracy_optimized'));
                 cd(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces','Accuracy_optimized'));
                 savevars(Decision_name,'Decision_surface',Decision_surface_acc_opt, 'num_paths',num_paths_acc_opt, 'Paths', Paths_acc_opt,'Spot_ID', Spot_ID_acc_opt); %will save the workspace into results.mat
             end 
              mkdir(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces','Accuracy'))
                 cd(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces','Accuracy'))
                 savevars(Decision_name,'Decision_surface',Decision_surface_acc, 'num_paths',num_paths_acc, 'Paths', Paths_acc,'Spot_ID', Spot_ID_acc); %will save the workspace into results.mat
                mkdir(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces','Performance')); 
                 cd(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces','Performance')); 
                 savevars(Decision_name,'Decision_surface',Decision_surface_perf,'num_paths',num_paths_perf,'Paths',Paths_perf, 'Spot_ID', Spot_ID_perf); %will save the workspace into results.mat
       else
             if HyperParameter_Search
              savevars(Decision_name,'Decision_surface',Decision_surface,'Decision_surface_optimized',Decision_surface_opt,'num_paths',num_paths_TF,'Paths', Paths_TF,'Spot_ID',Spot_ID_TF);
             else
              savevars(Decision_name,'Decision_surface',Decision_surface,'num_paths',num_paths_TF,'Paths', Paths_TF,'Spot_ID',Spot_ID_TF);
             end
       end    
     end
     cd(fullfile(CurrentFolder,'Classifiers','RCE','classify'))
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The END of the resamplings with k-fold each time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cd(fullfile(folder_name,'Rce_Results',discos))
if ~Estimate_size_time
  if HyperParameter_Search  
    [FinalConfmat,cum_perf_Final,acc_for_each_iteration,max_performance_per_cluster,performance_each_cluster_test,performance_each_cluster_train,...
    num_paths_for_max_clusters,num_paths_iter,num_cluster_iter,all_paths_freq,FinalConfmat_opt,cum_perf_Final_opt,acc_for_each_iteration_opt] = concat_variables(CurrentFolder,HyperParameter_Search, max_accuracy_Decision_Surface,no_classes,RCE_steps,kfold,topResample,paths,length(Hyperparameter));
    mean_acc = mean (acc_for_each_iteration);
    stdv_acc = std (acc_for_each_iteration);  
    mean_acc_opt = mean (acc_for_each_iteration_opt);
    stdv_acc_opt = std (acc_for_each_iteration_opt);
    cd(fullfile(folder_name,'Rce_Results',discos))
    save('acc_for_each_iteration_optim.mat','acc_for_each_iteration_opt')
    save('acc_for_each_iteration.mat','acc_for_each_iteration')
  else
    [FinalConfmat,cum_perf_Final,acc_for_each_iteration,max_performance_per_cluster,performance_each_cluster_test,performance_each_cluster_train,...
    num_paths_for_max_clusters,num_paths_iter,num_cluster_iter,all_paths_freq] = concat_variables(CurrentFolder,HyperParameter_Search, max_accuracy_Decision_Surface,no_classes,RCE_steps,kfold,topResample,paths);
    mean_acc = mean (acc_for_each_iteration);
    stdv_acc = std (acc_for_each_iteration); 
    cd(fullfile(folder_name,'Rce_Results',discos))
    save('acc_for_each_iteration.mat','acc_for_each_iteration')
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       replaced by acc_cal function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[accuracy,precision,recall,specificity]=acc_cal(cum_perf_Final,FinalConfmat);
precision=squeeze(precision);
recall=squeeze(recall);
specificity=squeeze(specificity);
num_paths=round(mean(num_paths_iter)); %The avarege
max_performance_each_cluster = mean(max_performance_per_cluster,1);
avg_performance_each_cluster_train = mean(performance_each_cluster_train,1);
avg_performance_each_cluster_test = mean(performance_each_cluster_test,1);
mean_num_paths_for_max_clusters = floor(mean(num_paths_for_max_clusters,1));
num_clusters = mean(num_cluster_iter,1);
res_out=[num_clusters;num_paths;accuracy;precision;recall;specificity]';

if HyperParameter_Search
[accuracy_opt,precision_opt,recall_opt,specificity_opt]=acc_cal(cum_perf_Final_opt,FinalConfmat_opt);
precision_opt=squeeze(precision_opt);
recall_opt=squeeze(recall_opt);
specificity_opt=squeeze(specificity_opt);
res_out_opt = [num_clusters;num_paths;accuracy_opt;precision_opt;recall_opt;specificity_opt]';
end
cd(fullfile(folder_name,'Rce_Results',discos))

fout =  fopen('Finalaccuracy_BasedOnScores.txt','w');
fprintf(fout,'#Clusters\t#Paths\tAcc\t');
for i=1:no_classes
fprintf(fout,'Precision_%s\tRecall_%s\tSpecificity_%s\t',classname{i},classname{i},classname{i});
end
for j=1:RCE_steps
    fprintf(fout,'\n%d\t%d\t%1.3%f',res_out(j,1),res_out(j,2),res_out(j,3));
   for k=1:no_classes
         fprintf(fout,'\t%1.3f\t%1.3f\t%1.3f',res_out(j,no_classes+3),res_out(j,2*no_classes+3),res_out(j,3*no_classes+3));
   end
end
fclose(fout);

if HyperParameter_Search
fout =  fopen('Finalaccuracy_BasedOnScores_paramoptim.txt','w');
fprintf(fout,'#Clusters\t#Paths\tAcc\t');
for i=1:no_classes
fprintf(fout,'Precision%s\tRecall%s\tSpecificity%s\t',classname{i},classname{i},classname{i});
end
for j=1:RCE_steps
    fprintf(fout,'\n%d\t%d\t%1.3%f',res_out_opt(j,1),res_out_opt(j,2),res_out_opt(j,3));
   for k=1:no_classes
         fprintf(fout,'\t%1.3f\t%1.3f\t%1.3f',res_out_opt(j,no_classes+3),res_out_opt(j,2*no_classes+3),res_out_opt(j,3*no_classes+3));
   end
end
fclose(fout);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Plot the track Accuracy
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% st={};  st_max={};
% for i=1:length(num_clusters)
%     str_cell{i} = mat2str(num_clusters(i));
%     str_paths{i} = mat2str(num_paths(i));
%     st{i} = strcat(mat2str(num_clusters(i)),'/',mat2str(num_paths(i)) );
%     st_max= mat2str(num_paths_for_max_clusters(i));
% end
% num = 1:1:length(num_clusters);
% for nc=1:no_classes
%     figure_1(nc,no_classes,classname(nc),accuracy,recall(nc,:),specificity(nc,:),precision(nc,:),...
%                  num_clusters,max_performance_each_cluster,avg_performance_each_cluster,...
%                  TOP_RESAMPLE,outer_kfold,st,num);
% end
% hold off;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Output the maximum accuracy over the clusters
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure_2(cum_perf_macro,counter,num_paths,accuracy,tp_rate,tn_rate,...
% 
% for nc=1:no_classes
%     figure_2(nc,classname(nc),cum_acc_macro,squeeze(cum_perf_macro(:,nc,:)),counter,num_paths,accuracy,recall(nc,:),specificity(nc,:),precision(nc,:),mean(acc_for_each_iteration),...
%           num_clusters,max_performance_each_cluster,...
%           avg_performance_each_cluster,st_max,TOP_RESAMPLE,resample,outer_kfold,st,num);
% end
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = sum(all_paths_freq,2);%sum each rows
[gsorted,gindex] = sort (freq,'descend');
significant_paths = paths(gindex);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save significant paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the location of the classes
total_freq = (topResample*kfold*length(num_paths));
n=length(groups_name);
grp=classname; 
no_classes = length(classname);
loc=zeros(no_classes,1);
for  j=1:no_classes
    for jj=1:n
        loc(j)=loc(j)+strcmp(groups_name(jj),grp{j});
    end
end                       
first=1;d = cell(no_classes,1);
for j=1:no_classes
    d{j}=data(first:first+loc(j)-1,:);
    first=first+loc(j);
end  

%fsamples = fopen('significant_paths_SVMKMEANS_samples.txt','w');
fsamples = fopen('significant_paths_KMEANS_samples.txt','w');
fid = fopen('significant_paths_KMEANS.txt','w');
fpaths = fopen('significant_paths_KMEANS_paths_levels.txt','w');
fprintf(fid,'SpotId\tPathName\t');
for i=1:no_classes
fprintf(fid,'%s\t',classname{i});
end
fprintf(fid,'Freq/%d\tScore',total_freq);

fprintf(fsamples,'SpotId\tpathName');
for i=1:n_samples
    fprintf(fsamples,'\t%s',mat2str(samples_name(i)));
end
fprintf(fsamples,'\tScore\n');
fprintf(fsamples,'SpotId\tpathName');
for i=1:n_samples
    fprintf(fsamples,'\t%s',mat2str(groups_ID(i)));
end
fprintf(fsamples,'\tScore');

fprintf(fpaths,'SpotId\tPathName');
for k =1:length(num_paths)
    st = sprintf('#paths:%d',num_paths(k));
    fprintf(fpaths,'\t%s',st);
end
sth = sprintf('\t%s/%d\t%s','Freq',total_freq,'Score');
fprintf(fpaths,sth);
for i=1:num_paths
    indx = gindex(i);
    score = (gsorted(i))/(topResample*kfold*length(num_paths));
    fprintf(fid,'\n%s\t%s',...
        mat2str(spot_id(indx)),mat2str(cell2mat(paths(indx))));
    for j=1:no_classes
        dat = d{j};
    fprintf(fid,'\t%1.3f',mean(dat(:,indx)));
    end
    clear dat
    fprintf(fid,'\t%d\t%2.3f',gsorted(i),score);
    fprintf(fsamples,'\n%s\t%s',...
        mat2str(spot_id(indx)),mat2str(cell2mat(paths(indx))));
    for kk=1:n_samples
        fprintf(fsamples,'\t%1.3f',data(kk,indx));
    end
    fprintf(fsamples,'\t%2.3f',score);

    %print for each paths its track at each level
    fprintf(fpaths,'\n%s\t%s',mat2str(spot_id(indx)),mat2str(cell2mat(paths(indx))) );
    for k =1:length(num_paths)
        fprintf(fpaths,'\t%d',all_paths_freq(indx,k) );
    end
    fprintf(fpaths,'\t%d\t%2.3f',gsorted(i),score );
end
fclose('all');
end

if Decision_Surface_Flag && Estimate_size_time 
if max_accuracy_Decision_Surface
   cd(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces','Accuracy'));
   S= dir('Decision_Surface1.mat');
   Decision_Surf_size = S.bytes;  
   cd(fullfile(CurrentFolder,'Rce_Results',discos,'Decision_Surfaces','Performance'));
   S= dir('Decision_Surface1.mat');
   Decision_Surf_size =  Decision_Surf_size+ S.bytes; 
else
   cd(fullfile(folder_name,'Rce_Results',discos,'Decision_Surfaces')); 
   S= dir('Decision_Surface1.mat');
   Decision_Surf_size = S.bytes; 
end
   varargout{2} = Decision_Surf_size;
end
if ~Estimate_size_time
   cd(fullfile(folder_name,'Rce_Results',discos))
    save('Results.mat', '-v7.3'); %will save the workspace into results.mat
end
time=toc(timer); varargout{1}=time;
 cd(fullfile(CurrentFolder,'Classifiers','RCE','classify'))
 clearvars Decision_surface* num_paths* Paths* Spot_ID* dat d i j k kk


