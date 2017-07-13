function [tp, fn, tn, fp, Conf_mat_bst,bestparam,varargout]= optimPara(data,group_no,data_test, group_no_test,test,kfold,classes,Decision_Surface_Flag,Hyperparameter)
% data: each arrow is one observation
% species : arrow(i) is the class member of abservation data(i)
% clsname : The name of the positive class or the case class

n_samples_tot= length(test);
samples_freq = zeros(1,n_samples_tot);% keep the frequency of the samples been tested
tp=0; tn=0; fp=0; fn=0;

cp2 = classperf(group_no_test);
itr= 5;
% Search for the best set of parameters
Optimizable_parameters =length(Hyperparameter);
Accuracy_paramset = zeros(prod(sizeout(Hyperparameter)),1);
Paramlist = paramgrid(Hyperparameter);
for i=1:prod(sizeout(Hyperparameter))
    Paramset =Paramlist(i,:);
    cp = classperf(group_no);
    for jj=1:itr
    indices = crossvalind('Kfold',group_no,kfold);
    for ij=1:kfold
         test = (indices == ij); train =~ test;
        Decision_Surf = classifier_train(data(train,:),group_no(train),Paramset);
        [classes_act] = classifier_predict(Decision_Surf, data(test,:));
        classperf(cp,classes_act,test);
    end
    end
    Accuracy_paramset(i)=cp.CorrectRate;
end
[max_Accuracy, Index] =  max(Accuracy_paramset);
bestparam= Paramlist(Index,:);

% Train using the best value of parameters
Decision_Surf_best = classifier_train(data, group_no, bestparam);
[classes_act_test] = classifier_predict(Decision_Surf_best,data_test);

% Updates the CP object with the current classification results
classperf(cp2,classes_act_test);

for j =1:length(test)
    if (test(j)==1)
        samples_freq(j)   = samples_freq(j)+1; %Keep the frequncy of the samples
    end
end

% Calculate performance Metrics
Conf_mat_bst=confusionmat(group_no_test,classes_act_test,'order',classes);
L=length(unique(classes));
for KJ=1:L
    tp(KJ)=Conf_mat_bst(KJ,KJ);
    fn(KJ)=sum(Conf_mat_bst(KJ,:))-Conf_mat_bst(KJ,KJ);
    fp(KJ)=sum(Conf_mat_bst(:,KJ))-Conf_mat_bst(KJ,KJ);
    tn(KJ)=sum(sum(Conf_mat_bst))-sum(Conf_mat_bst(KJ,:))-sum(Conf_mat_bst(:,KJ))+Conf_mat_bst(KJ,KJ);
end
if Decision_Surface_Flag == 1
    varargout{1}= Decision_Surf_best;
end 

function [no_para] = sizeout(x)
nout = length(x);no_para =zeros(1,nout);
for k = 1:nout
   no_para(k) = length(x{k});
end