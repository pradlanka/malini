function [max_avg_performance_each_cluster, Conf_mat_test, track_test_performance, num_clusters, itr_num_paths, paths_freq, varargout]= core_svmKmeansGeneSelection_v3(data,spot_id,paths,groups_name,samples_name,init_nmc,end_nmc,end_paths,decrease,kfold,resample,...
   classes,data_test,groups_test,samples_test,RCE_steps,Decision_Surface_Flag,HyperParameter_Search,varargin)
if nargin > 18
Hyperparameter = varargin{1};
end

fprintf('\nStart the Learning...');

DECFLAGE=0;         %decrease by percentage
if ( decrease>= 1)
    DECFLAGE = 1;   %decrease by integer number
end

%Initializing variables
max_avg_performance_each_cluster = zeros(4,RCE_steps);
Decision_surf = cell(1,RCE_steps);
track_test_performance=zeros(4,length(unique(classes)),RCE_steps);
Conf_mat_test = zeros(length(unique(classes)),length(unique(classes)),RCE_steps);

if HyperParameter_Search
Conf_mat_test_opt = zeros(length(unique(classes)),length(unique(classes)),RCE_steps);
track_test_performance_opt=zeros(4,length(unique(classes)),RCE_steps);
Decision_surf_opt = cell(1,RCE_steps);
end

n_paths = numel(paths);
nmc = init_nmc;  strt_nmc = nmc;
itr_num_paths=[];  nsrv_paths = size(data,1);
itr = 0;
num_clusters=[];        %keep the number of clusters at each iteration
num_paths=[];           %Keep track of the number of paths at each cluster
first_round=0; Spot_ID=[];
paths_idnt   = 1:1:n_paths;                               %Initialization
paths_levels = 1:1:n_paths;
paths_freq=[];  get_in = 0;  f_r = 1 ;   
while ( nmc >= end_nmc && nsrv_paths >= end_paths)
    first_round= first_round + 1;                           %Repeat the first round n times
    itr = itr +1 ;  num_clusters(itr)=nmc;
    [cidx, ctrs] = kmeans(data, nmc,'dist','correlation','rep',1,'disp','off','EmptyAction','singleton');
    results=zeros(1,nmc);
    track_test_performance_each_cluster=zeros(1,nmc);
    fprintf('\n');
    track_performance =zeros(nmc,RCE_steps);
    for c = 1:nmc
        n = sum(cidx==c);
        [track_performance(c,itr)]= runSVM_v3(data((cidx==c),:)',groups_name',kfold,10);
        prfc = track_performance(c,itr);
        [tp, fn, tn, fp, Conf_mat] = runSVMOnTestSet(data((cidx==c),:)',groups_name,data_test((cidx==c),:)',groups_test,samples_test,classes);
        track_test_performance_each_cluster(c) = sum(tp)/sum(sum(Conf_mat)); %the accuracy for each cluster based on the test part;
        results(c) = prfc; %track_test_performance_each_cluster(c);
    end
    %Train
    track_performance(strt_nmc+1,itr) = sum(track_performance(:,itr))/(nmc); %Average performance
    track_performance(strt_nmc+2,itr) = max(track_performance(:,itr));       %maximum Performance
    %Test
    [y indx] = sort(results); %Results from the TEST stage
    mc1 = indx(end);
    n_max = sum(cidx==mc1); %The number of paths in the cluster with maximum score based on the training
    max_acc = y(end);
    %Calculate the accuracy of the maximum cluster on the test samples
    max_avg_performance_each_cluster(1,itr)  = sum((track_test_performance_each_cluster))/nmc; %Average performance on the test data
    max_avg_performance_each_cluster(2,itr) = sum(results)/(nmc); % Average performance per cluster on the train dat
    max_avg_performance_each_cluster(3,itr)  = max_acc; %Train performance on the cluster with the maximum accurecy 
    max_avg_performance_each_cluster(4,itr ) = n_max ;  %the number of paths in this cluster

    [track_performance(strt_nmc+3,itr)] = runSVM_v3(data',groups_name',kfold,10);%Performance based on all the paths
    
    if ~Decision_Surface_Flag      
    [tp, fn, tn, fp, Conf_mat_perf] = runSVMOnTestSet(data',groups_name,data_test',groups_test,samples_test,classes,Decision_Surface_Flag); %Test Set
    else
    [tp, fn, tn, fp, Conf_mat_perf, Decision_surf{itr}] = runSVMOnTestSet(data',groups_name,data_test',groups_test,samples_test,classes,Decision_Surface_Flag); %Test Set
    end
        
    track_test_performance(1,:,itr) = tp ;    
    track_test_performance(2,:,itr) = fn ;
    track_test_performance(3,:,itr) = tn ;   
    track_test_performance(4,:,itr) = fp ;
    Conf_mat_test(:,:,itr)= Conf_mat_perf;
    
    if HyperParameter_Search
        if Decision_Surface_Flag 
            [tp_opt, fn_opt, tn_opt, fp_opt, Conf_mat_test_opt(:,:,itr), bestparam, Decision_surf_opt{itr}]= optimPara(data',groups_name,data_test',groups_test,samples_test,kfold,classes,Decision_Surface_Flag, Hyperparameter) ;
        else
            [tp_opt, fn_opt, tn_opt, fp_opt, Conf_mat_test_opt(:,:,itr),bestparam]= optimPara(data',groups_name,data_test',groups_test,samples_test,kfold,classes,Decision_Surface_Flag, Hyperparameter) ;
        end
    track_test_performance_opt(1,:,itr) = tp_opt ;    
    track_test_performance_opt(2,:,itr) = fn_opt ;
    track_test_performance_opt(3,:,itr) = tn_opt ;   
    track_test_performance_opt(4,:,itr) = fp_opt ;
    end

% Print the performance measures 
 if HyperParameter_Search
    fprintf('\n*Train:Clusters:%d #Paths:%d  Acc(all paths)=%f Avg=%f  Max=%f',nmc,nsrv_paths,track_performance(strt_nmc+3,itr),...
        track_performance(strt_nmc+1,itr),track_performance(strt_nmc+2,itr) );
    fprintf('\n#Test:Clusters:%d #Paths:%d Acc(all paths)=%f Avg=%f  Max(%d)=%f',nmc,nsrv_paths,sum(tp)/sum(sum(Conf_mat_perf)),...
        max_avg_performance_each_cluster(1,itr),n_max,max_avg_performance_each_cluster(2,itr)     );
    fprintf('\nTrue positives');disp(tp_opt); disp('    Confusionmatrix  ');disp(squeeze(Conf_mat_test_opt(:,:,itr)));
    
 else
       fprintf('\n*Train:Clusters:%d #Paths:%d  Acc(all paths)=%f Avg=%f  Max=%f',nmc,nsrv_paths,track_performance(strt_nmc+3,itr),...
        track_performance(strt_nmc+1,itr),track_performance(strt_nmc+2,itr) );
    fprintf('\n#Test:Clusters:%d #Paths:%d Acc(all paths)=%f Avg=%f  Max(%d)=%f',nmc,nsrv_paths,sum(tp)/sum(sum(Conf_mat_perf)),...
        max_avg_performance_each_cluster(1,itr),n_max,max_avg_performance_each_cluster(2,itr)     );
    fprintf('\nTrue positives');disp(tp); disp('    Confusionmatrix  ');disp(Conf_mat_perf);
 end
   
    track_performance(strt_nmc+4,itr)  = nsrv_paths;        %Number of paths
    itr_num_paths(itr) = nsrv_paths;                        %keep track about the number of paths at each iteration
    track_performance(strt_nmc+5,itr)  = nmc ;              %Number of clusters

    Paths{itr}=paths;
    num_paths(itr) = nsrv_paths;
    Spot_ID{itr}=spot_id;
    paths_freq(:,itr) = ismember(paths_idnt,paths_levels);

    [y indx] = sort(results);                               %sort the performance of each class
    if (DECFLAGE == 1)
        limit = decrease;
        nmc = nmc - decrease;
    else
        last_nmc = nmc;
        nmc = floor(nmc-decrease*nmc);
        limit = last_nmc - nmc;
    end;
    clust_index=[];
    for i=1:limit
        clust_index=[clust_index indx(i)];
    end
    if (nmc >=end_nmc )
        data=data(~ismember(cidx,clust_index),:);           %delete NOT significant paths
        data_test=data_test(~ismember(cidx,clust_index),:); %Remove them from TEST samples
        mask = ~ismember(cidx,clust_index);
        paths=paths(mask);
        spot_id=spot_id(mask);
        paths_levels=paths_levels(mask);
    end   
    nsrv_paths = size(data,1);
    if (f_r == 10)                                          %Only to avoid the repeated splits
        if (first_round<=15 )                               %To repeat the first round
            nmc = 10;                                       %Remember to change nmc up
        else
            f_r = 0 ;   get_in = 1;
        end
    end
    if (get_in == 1 )
        nmc = init_nmc; get_in = 0 ;
    end
end

if Decision_Surface_Flag  
varargout{1}= Decision_surf;
varargout{2} = num_paths; 
varargout{3} = Paths;
varargout{4} = Spot_ID;
if HyperParameter_Search 
varargout{5} = track_test_performance_opt;
varargout{6} = Conf_mat_test_opt;
varargout{7} = Decision_surf_opt;
varargout{8} = bestparam;
end
elseif HyperParameter_Search 
varargout{1} = track_test_performance_opt;
varargout{2} = Conf_mat_test_opt; 
varargout{3} = bestparam;
end
  